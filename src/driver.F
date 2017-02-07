!
      subroutine driver(nstep,istep)
!
! Driver layer of the model for current time step makes most calls
!   to dynamics and composition routines.
! Routines that do not have to do message-passing are generally called
!   from inside a latitude scan in this routine. Routines that do have
!   to make mpi calls are called once from this routine, and latitude
!   scan(s) are inside the called routine.
!
! -----------------------------------------------------------------------
!     driver.F called by advance.F, which also calls:
!     -- addiag.F. hdif1.F, hdif2.F
! -----------------------------------------------------------------------
! 8/16/06 btf: adapted for vtgcm from tiegcm dynamics.F.
! 8/16/06 btf: adapted for vtgcm from tiegcm dynamics.F.
!     ON:   compn2, cpktkm, swdot, chapman, qrj, hdif3, dt
!     OFF:  duv, comp_oco, comp,  
! -----------------------------------------------------------------------
! 11/21/06 swb: testing of sequential subroutines
!     ON:   compn2, cpktkm, swdot, chapman, qrj, hdif3, dt, duv, chemrates,
!           ionden, comp_oco, comp
!     OFF:  elden
! -----------------------------------------------------------------------
! 11/29/06 swb: testing of sequential subroutines
!     ON:   compn2, cpktkm, swdot, chapman, qrj, hdif3, dt, duv, chemrates,
!           ionden, comp_oco, comp, settei
!     OFF:  elden
! -----------------------------------------------------------------------
! 12/04/06 swb: testing of sequential subroutines
!     ON:   compn2, cpktkm, swdot, chapman, qrj, hdif3, dt, duv, chemrates,
!           ionden, comp_oco, comp, settei, eddyflds
!     OFF:  elden
! -----------------------------------------------------------------------
! 02/15/07 swb: testing of sequential subroutines
!     ON:   compn2, cpktkm, swdot, chapman, qrj, hdif3, dt, duv, chemrates,
!           ionden, comp_oco, comp, settei, eddyflds, comp_o2, minor_o2, minor
!     OFF:  elden
! -----------------------------------------------------------------------
! 03/12/08 swb: testing of sequential subroutines
! 04/16/08 swb: modify calls to comp_n4s and comp_n2d with n2p
!     ON:   compn2, cpktkm, swdot, chapman, qrj, hdif3, dt, duv, chemrates,
!           ionden, comp_oco, comp, settei, eddyflds, comp_o2, minor_o2, minor,
!           comp_n4S, minor_n4s, comp_n2d
!     OFF:  elden
! 05/12/08 swb: testing of sequential subroutines (comp_n2d, comp_no, comp_n4s)
!     ON:   compn2, cpktkm, swdot, chapman, qrj, hdif3, dt, duv, chemrates,
!           ionden, comp_oco, comp, settei, eddyflds, comp_o2, minor_o2, minor,
!           comp_n4S, minor_n4s, comp_n2d,  comp_no, minor_no
!     OFF:  elden
! 09/24/10 swb: testing of sequential subroutines (comp_so2, comp_so)
!     ON:   compn2, cpktkm, swdot, chapman, qrj, hdif3, dt, duv, chemrates,
!           ionden, comp_oco, comp, settei, eddyflds, comp_o2, minor_o2, minor,
!           comp_n4S, minor_n4s, comp_n2d,  comp_no, minor_no
!           --modified chapman.F and qrj.F for Jso, Jso2a,Jso2b,Jso3
!     OFF:  elden
! 06/06/14 swb: testing of aerosol heating in dt.F (add z argument to call stmt)
! 07/11/16 swb: testing of subroutines (comp_oh, chemrates, etc)
! -----------------------------------------------------------------------
!
      use fields_module
      use params_module,only: nlevp1
      use cons_module,only: kut
      use mpi_module,only: lon0,lon1,lat0,lat1
      use filter_module,only: filter_sp
      use addfld_module,only: addfld
      use hdif_module,only: hdif3,hdif_periodic
      use chapman_module,only: chapman
      use qrj_module,only: qrj
      use chemrates_module,only: chemrates_tdep,fs
      use dt_module,only: dt
      use ionden_module,only: ionden
      use eddyflds_module,only: eddyflds
      use o2_module,only: comp_o2, minor_o2, integral
      use n4s_module,only: comp_n4s, minor_n4s, comp_n2d,
     |                     comp_no, minor_no
      use so_module,only: comp_so, minor_so
      use so2_module,only: comp_so2, minor_so2
      implicit none
!
! Args:
      integer,intent(in) :: nstep,istep
!
! Local:
      integer :: i,j,k,n,lat,ier
      integer :: i0,i1,nk,nkm1,nlats,k0,k1
      logical,parameter :: debug=.false.   ! prints to stdout
      logical,parameter :: diags=.false.   ! calls to addfld
!
      if (debug) write(6,"(/,'Enter driver.')")
      i0=lon0 ; i1=lon1 ; nk=nlevp1 ; nkm1=nk-1 ! for addfsech
      k0=1 ; k1=nlevp1
      nlats = lat1-lat0+1
!
! This call moved to main tgcm.F
!     call bndcmp
!     if (debug) write(6,"('driver after bndcmp')")
!
! Latitude scan:
      do lat=lat0,lat1
!
! n2 composition:
        call compn2(
     |    tn  (levd0,lond0,lat,itp), ! 4d input
     |    barm(levd0,lond0,lat),     ! 3d input
     |    n2  (levd0,lond0,lat,itc), ! psi n2 output (mmr)
     |    1,nlevp1,lon0,lon1,lat)
!
!       if (diags) then
!         call addfld('N2_UNFIL',' ',' ',n2(:,lon0:lon1,lat,itc),
!    |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
!       endif
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('driver after compn2')")
!
! Filter n2:
      call filter_sp(n2(:,lon0:lon1,lat0:lat1,itc),
     |  1,nlevp1,lon0,lon1,lat0,lat1,kut)
      if (debug) write(6,"('driver after filter n2')")
!     if (diags) then
        do lat=lat0,lat1
          call addfld('N2_FILT',' ',' ',n2(:,lon0:lon1,lat,itc),
     |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
        enddo
!     endif
!     do lat=lat0,lat1
!       call addfld('Z_itc1',' ',' ',z(:,lon0:lon1,lat,itc),
!    |    'lev',1,nlevp1,'lon',lon0,lon1,lat)
!     enddo
!
! -------------------------------------------------------------------
! cp,kt,km:
! -------------------------------------------------------------------
!
      do lat=lat0,lat1
        call cpktkm(
     |    tn (levd0,lond0,lat,itp),
     |    o1 (levd0,lond0,lat,itp),
     |    co (levd0,lond0,lat,itp),
     |    co2(levd0,lond0,lat),
     |    n2 (levd0,lond0,lat,itc),
     |    cp (levd0,lond0,lat),      ! output
     |    kt (levd0,lond0,lat),      ! output
     |    km (levd0,lond0,lat),      ! output
     |    1,nlevp1,lon0,lon1,lat)

        if (diags) then
          call addfld('CP',' ',' ',cp(:,lon0:lon1,lat),
     |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
          call addfld('KT',' ',' ',kt(:,lon0:lon1,lat),
     |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
          call addfld('KM',' ',' ',km(:,lon0:lon1,lat),
     |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
        endif
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('driver after cpktkm')")
!     do lat=lat0,lat1
!       call addfld('Z_itc2',' ',' ',z(:,lon0:lon1,lat,itc),
!    |    'lev',1,nlevp1,'lon',lon0,lon1,lat)
!     enddo
!
! -----------------------------------------------------------------
! Calculate omega for vertical velocity:
! (latitude scan is inside swdot)
! -----------------------------------------------------------------
!
      if (debug) write(6,"('driver before swdot')")
      call swdot(
     |  un(levd0,lond0,latd0,itp), ! un input
     |  vc(levd0,lond0,latd0,itp), ! vc input (from addiag)
     |  w (levd0,lond0,latd0,itc), ! omega output
     |  1,nlevp1,lon0,lon1,lat0,lat1,lat0,lat1)
      if (debug) write(6,"('driver after swdot')")
      if (diags) then
        do lat=lat0,lat1
          call addfld('W_UPD',' ',' ',w(:,lon0:lon1,lat,itc),
     |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
        enddo
      endif
!     do lat=lat0,lat1
!       call addfld('Z_itc3',' ',' ',z(:,lon0:lon1,lat,itc),
!    |    'lev',1,nlevp1,'lon',lon0,lon1,lat)
!     enddo
!
! -----------------------------------------------------------------
! Calculate column densities and line integrals of o1, co, co2, n2
! Prognostics z,tn,o1,co,barm are input, diagnostics vo1,vco,vco2,
!   vn2,sco1,scco,scco2,scn2 are output.
! 31-args (standard suite for so,so_nm,so2,so2_nm on history)
! 33-args (start suite for so,so2 only on history; not so_nm, so2_nm)
! -----------------------------------------------------------------
!
      do lat=lat0,lat1
        call chapman(
     |    z    (levd0,lond0,lat,itc), ! updated Z from addiag
     |    tn   (levd0,lond0,lat,itp), ! 4d in
     |    o1   (levd0,lond0,lat,itp), ! 4d in
     |    co   (levd0,lond0,lat,itp), ! 4d in
     |    co2  (levd0,lond0,lat),     ! is local to chapman
     |    n2   (levd0,lond0,lat,itp), ! 4d in
     |    o2   (levd0,lond0,lat,itp), ! 4d in
     |    so   (levd0,lond0,lat,itp), ! 4d inout
     |    so2  (levd0,lond0,lat,itp), ! 4d inout
!    |    so_nm  (levd0,lond0,lat,itp), ! 4d inout
!    |    so2_nm (levd0,lond0,lat,itp), ! 4d inout
     |    barm (levd0,lond0,lat),     ! 3d in
     |    vo1  (levd0,lond0,lat),     ! 3d out
     |    vco  (levd0,lond0,lat),     ! 3d out
     |    vco2 (levd0,lond0,lat),     ! 3d out
     |    vn2  (levd0,lond0,lat),     ! 3d out
     |    vo2  (levd0,lond0,lat),     ! 3d out
     |    vso  (levd0,lond0,lat),     ! 3d out
     |    vso2 (levd0,lond0,lat),     ! 3d out
     |    vso3 (levd0,lond0,lat),     ! 3d out
     |    sco1 (levd0,lond0,lat),     ! 3d out
     |    scco (levd0,lond0,lat),     ! 3d out
     |    scco2(levd0,lond0,lat),     ! 3d out
     |    scn2 (levd0,lond0,lat),     ! 3d out
     |    sco2 (levd0,lond0,lat),     ! 3d out
     |    scso (levd0,lond0,lat),     ! 3d out
     |    scso2(levd0,lond0,lat),     ! 3d out
     |    scso3(levd0,lond0,lat),     ! 3d out
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('driver after chapman')")
!
! Post-processors cannot read 1.e80, so change these to
! 1.e36 if saving on secondary histories (then change back
! for rest of model):
      if (diags) then
        do lat=lat0,lat1
          do k=k0,k1
            do i=i0,i1
              if (sco2(k,i,lat)==1.e80) sco2(k,i,lat)=1.e36
              if (sco1(k,i,lat)==1.e80) sco1(k,i,lat)=1.e36
              if (scn2(k,i,lat)==1.e80) scn2(k,i,lat)=1.e36
              if (scso(k,i,lat)==1.e80) scso(k,i,lat)=1.e36
              if (scso2(k,i,lat)==1.e80) scso2(k,i,lat)=1.e36
              if (scso3(k,i,lat)==1.e80) scso3(k,i,lat)=1.e36
            enddo
          enddo
          call addfld('SCO1',' ',' ',sco1(:,lon0:lon1,lat),
     |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
          call addfld('SCCO',' ',' ',scco(:,lon0:lon1,lat),
     |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
          call addfld('SCCO2',' ',' ',scco2(:,lon0:lon1,lat),
     |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
          call addfld('SCN2',' ',' ',scn2(:,lon0:lon1,lat),
     |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
          call addfld('SCO2',' ',' ',sco2(:,lon0:lon1,lat),
     |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
          call addfld('SCSO',' ',' ',scso(:,lon0:lon1,lat),
     |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
          call addfld('SCSO2',' ',' ',scso2(:,lon0:lon1,lat),
     |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
          call addfld('SCSO3',' ',' ',scso3(:,lon0:lon1,lat),
     |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
          do k=k0,k1
            do i=i0,i1
              if (sco2(k,i,lat)==1.e36) sco2(k,i,lat)=1.e80
              if (sco1(k,i,lat)==1.e36) sco1(k,i,lat)=1.e80
              if (scn2(k,i,lat)==1.e36) scn2(k,i,lat)=1.e80
              if (scso(k,i,lat)==1.e36) scso(k,i,lat)=1.e80
              if (scso2(k,i,lat)==1.e36) scso2(k,i,lat)=1.e80
              if (scso3(k,i,lat)==1.e36) scso3(k,i,lat)=1.e80
            enddo
          enddo
        enddo ! lat=lat0,lat1
      endif
!     do lat=lat0,lat1
!       call addfld('Z_itc4',' ',' ',z(:,lon0:lon1,lat,itc),
!    |    'lev',1,nlevp1,'lon',lon0,lon1,lat)
!     enddo
!
! -----------------------------------------------------------------
! Calculate Te and Ti (from Fox and Sung (2001)) 
! Te,Ti needed in chemrates (only)
!     subroutine settei(tn,z,te,ti,lev0,lev1,lon0,lon1,lat)
! -----------------------------------------------------------------
!
      do lat=lat0,lat1
        call settei(
     |    tn  (levd0,lond0,lat,itp),
     |    z   (levd0,lond0,lat,itc), ! updated Z from addiag
     |    te  (levd0,lond0,lat),     ! output (new 3-D field)
     |    ti  (levd0,lond0,lat),     ! output (new 3-D field)
     |    1,nlevp1,lon0,lon1,lat)
      enddo
      if (debug) write(6,"('driver after settei')")
!
! -----------------------------------------------------------------
! Calculate temperature dependent reaction rates (chemrates_module):
! (Pass te,ti from settei)
! -----------------------------------------------------------------
!
      do lat=lat0,lat1
        call chemrates_tdep(
     |    tn  (levd0,lond0,lat,itp),
     |    te  (levd0,lond0,lat),
     |    ti  (levd0,lond0,lat),
     |    1,nlevp1,lon0,lon1,lat)
      enddo
      if (debug) write(6,"('driver after chemrates_tdep')")
!
! -----------------------------------------------------------------
! Calculate qrj fields: QEUV, JCO2, JO2, JN2, IO, ICO, IN2, ICO2, IO2
!                             JSO, JSO2a, JSO2b, JSO3
! -----------------------------------------------------------------
      do lat=lat0,lat1
        call qrj(
     |    sco1(levd0,lond0,lat),
     |    scco(levd0,lond0,lat),
     |    scn2(levd0,lond0,lat),
     |    sco2(levd0,lond0,lat),
     |    scco2(levd0,lond0,lat),
     |    scso(levd0,lond0,lat),
     |    scso2(levd0,lond0,lat),
     |    scso3(levd0,lond0,lat),
     |    tn(levd0,lond0,lat,itp),
     |    o1(levd0,lond0,lat,itp),
     |    co(levd0,lond0,lat,itp),
     |    n2(levd0,lond0,lat,itp),
     |    o2(levd0,lond0,lat,itp),
     |    so(levd0,lond0,lat,itp), 
     |    so2(levd0,lond0,lat,itp),
     |    barm(levd0,lond0,lat),
     |    cp(levd0,lond0,lat),
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1 for qrj
      if (debug) write(6,"('driver after qrj')")
!
! -----------------------------------------------------------------
! Make hdt,u,v,o1,co:
! -----------------------------------------------------------------
!
      do lat=lat0,lat1
        call hdif3(
     |    cp(levd0,lond0,lat),     ! specific heat input
     |    kldt ,                   ! 3d input from hdif2 for tn
     |    kldu ,                   ! 3d input from hdif2 for un
     |    kldv ,                   ! 3d input from hdif2 for vn
     |    kldo1,                   ! 3d input from hdif2 for o1
     |    kldco,                   ! 3d input from hdif2 for co
     |    hdt (levd0,lond0,lat),   ! 2d tn output at current lat
     |    hdu (levd0,lond0,lat),   ! 2d un output at current lat
     |    hdv (levd0,lond0,lat),   ! 2d vn output at current lat
     |    hdo1(levd0,lond0,lat),   ! 2d o1 output at current lat
     |    hdco(levd0,lond0,lat),   ! 2d h2 output at current lat
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('driver after hdif3')")
!
! -----------------------------------------------------------------
! Calculate ions using the "old" vtgcm2 method (prod/loss):
! This will be a diagnostic routine, i.e., to verify elden.
! Solve for electron density and ion species iteratively:
!
! -----------------------------------------------------------------
      do lat=lat0,lat1
        call ionden(
     |    tn  (levd0,lond0,lat,itp),
     |    barm(levd0,lond0,lat),
     |    o1  (levd0,lond0,lat,itp),
     |    co  (levd0,lond0,lat,itp),
     |    co2 (levd0,lond0,lat),
     |    n4s (levd0,lond0,lat,itp),
     |    n2  (levd0,lond0,lat,itp),
     |    n2p (levd0,lond0,lat),     ! output
     |    nop (levd0,lond0,lat),     ! output
     |    co2p(levd0,lond0,lat),     ! output
     |    o2p (levd0,lond0,lat,itc), ! output
     |    op  (levd0,lond0,lat,itc), ! output
     |    ne  (levd0,lond0,lat,itc), ! output
     |    1,nlevp1,lon0,lon1,lat)
      enddo
      if (debug) write(6,"('driver after ionden')")
!
! Solve for electron density and ion species using quart solver:
!     subroutine elden(tn,barm,o1,co,co2,n4s,n2,
!    |  n2p,nop,co2p,o2p,electrons,lev0,lev1,lon0,lon1,lat)
!
!     do lat=lat0,lat1
!       call elden(
!    |    tn  (levd0,lond0,lat,itp),
!    |    barm(levd0,lond0,lat),
!    |    o1  (levd0,lond0,lat,itp),
!    |    co  (levd0,lond0,lat,itp),
!    |    co2 (levd0,lond0,lat),
!    |    n4s (levd0,lond0,lat,itp),
!    |    n2  (levd0,lond0,lat,itp),
!    |    n2p (levd0,lond0,lat),     ! output
!    |    nop (levd0,lond0,lat),     ! output
!    |    co2p(levd0,lond0,lat),     ! output
!    |    o2p (levd0,lond0,lat,itc), ! output
!    |    op  (levd0,lond0,lat,itc), ! output
!    |    ne  (levd0,lond0,lat,itc), ! output
!    |    1,nlevp1,lon0,lon1,lat)
!     enddo
!
! -----------------------------------------------------------------
! Sources and sinks for major species composition O, CO:
!     subroutine comp_oco(tn,barm,o1,co,co2,n2,o2,op,co2p,o2p,ne,
!    |  lev0,lev1,lon0,lon1,lat)
! -----------------------------------------------------------------
      do lat=lat0,lat1
        if (debug) write(6,"('driver calling comp_oco: lat=',i3)") lat
        call comp_oco(
     |    tn   (levd0,lond0,lat,itp), ! neutral temperature
     |    barm (levd0,lond0,lat),     ! mean molecular weight
     |    o1   (levd0,lond0,lat,itp), ! atomic oxygen (mmr)
     |    co   (levd0,lond0,lat,itp), ! carbon momoxide (mmr)
     |    co2  (levd0,lond0,lat),     ! carbon dioxide (mmr)
     |    n2   (levd0,lond0,lat,itc), ! updated psi n2 from compn2
     |    o2   (levd0,lond0,lat,itp), ! molecular oxygen (mmr)
     |    op   (levd0,lond0,lat,itc), ! O+ (on history)
     |    co2p (levd0,lond0,lat),     ! CO2+ (diagnostic)
     |    o2p  (levd0,lond0,lat,itc), ! O2+ (on history)
     |    ne   (levd0,lond0,lat,itc), ! electron density (on history)
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('driver after comp_oco')")
!
! --------------------------------------------------------------------------
! Calculate eddy fields required for comp, dt, duv, minor
!     subroutine eddyflds(tn,barm,difk,dift,xmue,lev0,lev1,lon0,lon1,lat)
! -----------------------------------------------------------------
!
      do lat=lat0,lat1
        call eddyflds(
     |    tn  (levd0,lond0,lat,itp),  ! neutral temperature
     |    barm (levd0,lond0,lat),     ! mean molecular weight
     |    difk (levd0,lond0,lat),     ! Eddy Diffusion Coefficient (1/sec)
     |    dift (levd0,lond0,lat),     ! Eddy Thermal Coefficient (1/sec)
     |    xmue (levd0,lond0,lat),     ! Eddy Viscosity Coefficient (1/sec)
     |    1,nlevp1,lon0,lon1,lat)
      enddo
      if (debug) write(6,"('driver after eddyflds')")
!
! --------------------------------------------------------------------------
!     subroutine comp(tn,difk,o1,o1_nm,co,co_nm,n2,un,vn,w,hdo1,hdco,
!    |                o1_upd,o1nm_upd,co_upd,conm_upd,
!    |                lev0,lev1,lon0,lon1,lat0,lat1)
!
! Update main composition species O, CO: On (11/21/06)
      call comp(
     |  tn    (levd0,lond0,latd0,itp), ! neutral temperature
     |  difk  (levd0,lond0,latd0),     ! eddy diffusion coefficient (eddyflds)
     |  o1    (levd0,lond0,latd0,itp), ! O (mmr)
     |  o1_nm (levd0,lond0,latd0,itp), ! O (mmr) at time n-1
     |  co    (levd0,lond0,latd0,itp), ! CO (mmr)
     |  co_nm (levd0,lond0,latd0,itp), ! CO (mmr) at time n-1
     |  n2    (levd0,lond0,latd0,itc), ! updated n2 from compn2
     |  un    (levd0,lond0,latd0,itp), ! zonal velocity
     |  vn    (levd0,lond0,latd0,itp), ! meridional velocity
     |  w     (levd0,lond0,latd0,itp), ! vertical velocity
     |  hdo1  (levd0,lond0,latd0),     ! O horizontal diffusion (hdif3)
     |  hdco  (levd0,lond0,latd0),     ! CO horizontal diffusion (hdif3)
     |  o1    (levd0,lond0,latd0,itc), ! output: O updated for current step
     |  o1_nm (levd0,lond0,latd0,itc), ! output: O updated for previous step
     |  co    (levd0,lond0,latd0,itc), ! output: CO updated for current step
     |  co_nm (levd0,lond0,latd0,itc), ! output: CO updated for previous step
     |  1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('driver after comp')")
!     if (diags) then
!       do lat=lat0,lat1
!         call addfld('O1_UPD',' ',' ',o1(:,lon0:lon1,lat,itc),
!    |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
!         call addfld('CO_UPD',' ',' ',co(:,lon0:lon1,lat,itc),
!    |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
!       enddo
!     endif
!
! --------------------------------------------------------------------------
! Periodic points for horizontal diffusion terms output from hdif3:
! This may not be necessary.
      call hdif_periodic(hdt,hdu,hdv,hdo1,hdco,
     |    1,nlevp1,lon0,lon1,lat0,lat1)
! --------------------------------------------------------------------------
!
! Advance TN: ARG list modified on 6/06/14: swb
      if (debug) write(6,"('driver calling dt')")
      call dt(
     |  tn    (levd0,lond0,latd0,itp),
     |  tn_nm (levd0,lond0,latd0,itp),
     |  un    (levd0,lond0,latd0,itp),
     |  vn    (levd0,lond0,latd0,itp),
     |  o1    (levd0,lond0,latd0,itp),
     |  co    (levd0,lond0,latd0,itp),
     |  barm  (levd0,lond0,latd0),      ! mean molecular weight
     |  dift  (levd0,lond0,latd0),      ! eddy thermal conductivity (eddyflds)
     |  cp    (levd0,lond0,latd0),      ! specific heat        (cpktkm.F)
     |  kt    (levd0,lond0,latd0),      ! molecular diffusion  (cpktkm.F)
     |  km    (levd0,lond0,latd0),      ! molecular viscosity  (cpktkm.F)
     |  hdt   (levd0,lond0,latd0),      ! horizontal diffusion (hdif.F)
     |  w     (levd0,lond0,latd0,itc),  ! updated W (swdot.F)
     |  tn    (levd0,lond0,latd0,itc),  ! output updated tn
     |  tn_nm (levd0,lond0,latd0,itc),  ! output updated tn at time n-1
     |  1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('driver after dt')")
!     if (diags) then
!       do lat=lat0,lat1
!         call addfld('TN_UPD',' ',' ',tn(:,lon0:lon1,lat,itc),
!    |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
!         call addfld('TNNM_UPD',' ',' ',tn_nm(:,lon0:lon1,lat,itc),
!    |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
!       enddo
!     endif
!
! --------------------------------------------------------------------------
! Advance U,V:
      call duv(
     |  tn    (levd0,lond0,latd0,itp), ! neutral temperature
     |  tn    (levd0,lond0,latd0,itc), ! updated neutral temperature (dt.F)
     |  tn_nm (levd0,lond0,latd0,itp), ! tn at time n-1
     |  un    (levd0,lond0,latd0,itp), ! zonal velocity
     |  vn    (levd0,lond0,latd0,itp), ! meridional velocity
     |  un_nm (levd0,lond0,latd0,itp), ! zonal velocity at time n-1
     |  vn_nm (levd0,lond0,latd0,itp), ! meridional velocity at time n-1
     |  w     (levd0,lond0,latd0,itc), ! updated vertical velocity (swdot.F)
     |  barm  (levd0,lond0,latd0),     ! mean molecular weight
     |  xmue  (levd0,lond0,latd0),     ! eddy viscosity conductivity (eddyflds)
     |  z     (levd0,lond0,latd0,itp), ! geopotential height
     |  hdu   (levd0,lond0,latd0),     ! horizontal diffusion of U (hdif3.F)
     |  hdv   (levd0,lond0,latd0),     ! horizontal diffusion of V (hdif3.F)
!    |  ui    (levd0,lond0,latd0),     ! zonal ion velocity
!    |  vi    (levd0,lond0,latd0),     ! meridional ion velocity
!    |  lxx   (levd0,lond0,latd0),     ! xx ion drag coefficient
!    |  lyy   (levd0,lond0,latd0),     ! yy ion drag coefficient
!    |  lxy   (levd0,lond0,latd0),     ! xy ion drag coefficient
!    |  lyx   (levd0,lond0,latd0),     ! yx ion drag coefficient
     |  km    (levd0,lond0,latd0),     ! molecular viscosity  (cpktkm.F)
     |  un    (levd0,lond0,latd0,itc), ! output updated un
     |  un_nm (levd0,lond0,latd0,itc), ! output updated un at time n-1
     |  vn    (levd0,lond0,latd0,itc), ! output updated vn
     |  vn_nm (levd0,lond0,latd0,itc), ! output updated vn at time n-1
     |  1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('driver after duv')")
!     if (diags) then
!       do lat=lat0,lat1
!         call addfld('UN_UPD',' ',' ',un(:,lon0:lon1,lat,itc),
!    |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
!         call addfld('VN_UPD',' ',' ',vn(:,lon0:lon1,lat,itc),
!    |      'lev',1,nlevp1,'lon',lon0,lon1,lat)
!       enddo
!     endif
! --------------------------------------------------------------------------
! Advance o2 (minor):
!     subroutine comp_o2(tn,o1,co,co2,n2,o2,barm,o2p,op,co2p,ne,
!    |  lev0,lev1,lon0,lon1,lat)
! -----------------------------------------------------------------
      do lat=lat0,lat1
        if (debug) write(6,"('driver calling compo2: lat=',i3)") lat
        call comp_o2(
     |    tn   (levd0,lond0,lat,itc), ! neutral temperature
     |    o1   (levd0,lond0,lat,itc), ! atomic oxygen (mmr)
     |    co   (levd0,lond0,lat,itc), ! carbon momoxide (mmr)
     |    co2  (levd0,lond0,lat),     ! carbon dioxide (mmr)
     |    n2   (levd0,lond0,lat,itc), ! updated n2 (mmr) from compn2
     |    o2   (levd0,lond0,lat,itp), ! O2 (old, on history)
     |    barm (levd0,lond0,lat),     ! mean molecular weight
     |    o2p  (levd0,lond0,lat,itc), ! O2+ (on history)
     |    op   (levd0,lond0,lat,itc), ! O+ (on history)
     |    co2p (levd0,lond0,lat),     ! CO2+ (diagnostic)
     |    ne   (levd0,lond0,lat,itc), ! electron density (on history)
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('driver after comp_o2')")
! --------------------------------------------------------------------------
!     subroutine minor_o2(tn,difk,o1,co,co2,n2,o2,o2_nm1,o2_out,o2_nm1_out,
!    |  lev0,lev1,lon0,lon1,lat0,lat1)
! --------------------------------------------------------------------------
      call minor_o2(
     |  tn    (levd0,lond0,latd0,itc),   ! neutral temperature
     |  difk  (levd0,lond0,latd0),       ! eddy diffusion coefficient (eddyflds)
     |  o1    (levd0,lond0,latd0,itc),   ! atomic oxygen (mmr)
     |  co    (levd0,lond0,latd0,itc),   ! carbon momoxide (mmr)
     |  co2   (levd0,lond0,latd0),       ! carbon dioxide (mmr)
     |  n2    (levd0,lond0,latd0,itc),   ! updated n2 (mmr) from compn2
     |  o2    (levd0,lond0,latd0,itp),   ! previous O2 (at n) on history
     |  o2_nm (levd0,lond0,latd0,itp),   ! previous O2 (at n-1) on history
     |  o2    (levd0,lond0,latd0,itc),   ! new O2 (at n) on history
     |  o2_nm (levd0,lond0,latd0,itc),   ! new O2 (at n-1) on history
     |  1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('driver after minor_o2')")
! -----------------------------------------------------------------
! Advance n2d (pce species, no minor called): 15 inputs call
!     subroutine comp_n2d(tn,o1,co,co2,n2,barm,nop,n2p,ne,n2d,
!    |  lev0,lev1,lon0,lon1,lat)
! -----------------------------------------------------------------
      do lat=lat0,lat1
        if (debug) write(6,"('driver calling compn2d: lat=',i3)") lat
        call comp_n2d(
     |    tn   (levd0,lond0,lat,itc), ! neutral temperature
     |    o1   (levd0,lond0,lat,itc), ! atomic oxygen (mmr)
     |    co   (levd0,lond0,lat,itc), ! carbon momoxide (mmr)
     |    co2  (levd0,lond0,lat),     ! carbon dioxide (mmr)
     |    n2   (levd0,lond0,lat,itc), ! updated psi n2 from compn2
     |    barm (levd0,lond0,lat),     ! mean molecular weight
     |    nop  (levd0,lond0,lat),     ! NO+ (diagnostic)
     |    n2p  (levd0,lond0,lat),     ! N2+ (diagnostic)
     |    ne   (levd0,lond0,lat,itc), ! electron density (on history)
     |    n2d  (levd0,lond0,lat,itc), ! n2d output (mmr)
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('driver after comp_n2d')")
! --------------------------------------------------------------------------
! Advance no (minor): 15-inputs to call
!     subroutine comp_no(tn,o1,co,co2,n2,barm,n4s,n2d,no,o2p,
!    |  lev0,lev1,lon0,lon1,lat)
! -----------------------------------------------------------------
      do lat=lat0,lat1
        if (debug) write(6,"('driver calling comp_n4s: lat=',i3)") lat
        call comp_no(
     |    tn   (levd0,lond0,lat,itc), ! neutral temperature
     |    o1   (levd0,lond0,lat,itc), ! atomic oxygen (mmr)
     |    co   (levd0,lond0,lat,itc), ! carbon momoxide (mmr)
     |    co2  (levd0,lond0,lat),     ! carbon dioxide (mmr)
     |    n2   (levd0,lond0,lat,itc), ! updated psi n2 from compn2
     |    barm (levd0,lond0,lat),     ! mean molecular weight
     |    n4s  (levd0,lond0,lat,itp), ! n4s (old, on history)
     |    n2d  (levd0,lond0,lat,itc), ! n2d (current, from comp_n2d, on history)
     |    no   (levd0,lond0,lat,itp), ! no (old, on history)
     |    o2p  (levd0,lond0,lat,itc), ! O2+ (on history)
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('driver after comp_n4s')")
! --------------------------------------------------------------------------
!     subroutine minor_no(tn,difk,o1,co,co2,n2,no,no_nm1,no_out,no_nm1_out,
!    |  lev0,lev1,lon0,lon1,lat0,lat1)
! --------------------------------------------------------------------------
        call minor_no(
     |    tn   (levd0,lond0,latd0,itc),      ! neutral temperature
     |    difk (levd0,lond0,latd0),          ! eddy diffusion coefficient (eddyflds)
     |    o1   (levd0,lond0,latd0,itc),      ! atomic oxygen (mmr)
     |    co   (levd0,lond0,latd0,itc),      ! carbon momoxide (mmr)
     |    co2  (levd0,lond0,latd0),          ! carbon dioxide (mmr)
     |    n2   (levd0,lond0,latd0,itc),      ! updated n2 (mmr) from compn2
     |    no   (levd0,lond0,latd0,itp),      ! previous NO (at n) on history
     |    no_nm (levd0,lond0,latd0,itp),     ! previous NO (at n-1) on history
     |    no   (levd0,lond0,latd0,itc),      ! new NO (at n) on history
     |    no_nm (levd0,lond0,latd0,itc),     ! new NO (at n-1) on history
     |    1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('driver after minor_no')")
! --------------------------------------------------------------------------
! Advance n4s (minor): 19-inputs to call
!     subroutine comp_n4s(tn,o1,co,co2,n2,barm,n4s,n2d,no,o2p,op,nop,n2p,ne,
!    |  lev0,lev1,lon0,lon1,lat)
! -----------------------------------------------------------------
      do lat=lat0,lat1
        if (debug) write(6,"('driver calling comp_n4s: lat=',i3)") lat
        call comp_n4s(
     |    tn   (levd0,lond0,lat,itc), ! neutral temperature
     |    o1   (levd0,lond0,lat,itc), ! atomic oxygen (mmr)
     |    co   (levd0,lond0,lat,itc), ! carbon momoxide (mmr)
     |    co2  (levd0,lond0,lat),     ! carbon dioxide (mmr)
     |    n2   (levd0,lond0,lat,itc), ! updated psi n2 from compn2
     |    barm (levd0,lond0,lat),     ! mean molecular weight
     |    n4s  (levd0,lond0,lat,itp), ! n4s (old, on history)
     |    n2d  (levd0,lond0,lat,itc), ! n2d (current, from comp_n2d, on history)
     |    no   (levd0,lond0,lat,itc), ! no  (current, from comp_no,  on history)
     |    o2p  (levd0,lond0,lat,itc), ! O2+ (on history)
     |    op   (levd0,lond0,lat,itc), ! O+ (on history)
     |    nop  (levd0,lond0,lat),     ! NO+ (diagnostic)
     |    n2p  (levd0,lond0,lat),     ! N2+ (diagnostic)
     |    ne   (levd0,lond0,lat,itc), ! electron density (on history)
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('driver after comp_n4s')")
! --------------------------------------------------------------------------
!     subroutine minor_n4s(tn,difk,o1,co,co2,n2,n4s,n4s_nm1,n4s_out,n4s_nm1_out,
!    |  lev0,lev1,lon0,lon1,lat0,lat1)
! --------------------------------------------------------------------------
        call minor_n4s(
     |    tn   (levd0,lond0,latd0,itc),      ! neutral temperature
     |    difk (levd0,lond0,latd0),          ! eddy diffusion coefficient (eddyflds)
     |    o1   (levd0,lond0,latd0,itc),      ! atomic oxygen (mmr)
     |    co   (levd0,lond0,latd0,itc),      ! carbon momoxide (mmr)
     |    co2  (levd0,lond0,latd0),          ! carbon dioxide (mmr)
     |    n2   (levd0,lond0,latd0,itc),      ! updated n2 (mmr) from compn2
     |    n4s  (levd0,lond0,latd0,itp),      ! previous N4S (at n) on history
     |    n4s_nm (levd0,lond0,latd0,itp),    ! previous N4S (at n-1) on history
     |    n4s   (levd0,lond0,latd0,itc),     ! new N4S (at n) on history
     |    n4s_nm (levd0,lond0,latd0,itc),    ! new N4S (at n-1) on history
     |    1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('driver after minor_n4s')")
! -----------------------------------------------------------------
! Must populate so2_nm and so_nm since empty from starting history : test
      so2_nm(:,:,:,itp) = so2(:,:,:,itp)
      so_nm (:,:,:,itp) = so(:,:,:,itp)
! --------------------------------------------------------------------------
! Advance so (minor): 14-inputs to call
!     subroutine comp_so(tn,o1,co,co2,n2,o2,so,so2,barm,
!    |  lev0,lev1,lon0,lon1,lat)
! -----------------------------------------------------------------
      do lat=lat0,lat1
        if (debug) write(6,"('driver calling comp_so: lat=',i3)") lat
        call comp_so(
     |    tn   (levd0,lond0,lat,itc), ! neutral temperature
     |    o1   (levd0,lond0,lat,itc), ! atomic oxygen (mmr)
     |    co   (levd0,lond0,lat,itc), ! carbon momoxide (mmr)
     |    co2  (levd0,lond0,lat),     ! carbon dioxide (mmr)
     |    n2   (levd0,lond0,lat,itc), ! updated psi n2 from compn2
     |    o2   (levd0,lond0,lat,itc), ! updated psi o2 from compo2
     |    so   (levd0,lond0,lat,itp), ! so (old, on history)
!    |    so2  (levd0,lond0,lat,itc), ! so2 (current, from comp_so2, on history)
     |    so2  (levd0,lond0,lat,itp), ! so2 (old, on history)
     |    barm (levd0,lond0,lat),     ! mean molecular weight
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('driver after comp_so')")
! --------------------------------------------------------------------------
!     subroutine minor_so(tn,difk,o1,co,co2,n2,so,so_nm1,so_out,so_nm1_out,
!    |  lev0,lev1,lon0,lon1,lat0,lat1)
! --------------------------------------------------------------------------
        call minor_so(
     |    tn   (levd0,lond0,latd0,itc),      ! neutral temperature
     |    difk (levd0,lond0,latd0),          ! eddy diffusion coefficient (eddyflds)
     |    o1   (levd0,lond0,latd0,itc),      ! atomic oxygen (mmr)
     |    co   (levd0,lond0,latd0,itc),      ! carbon momoxide (mmr)
     |    co2  (levd0,lond0,latd0),          ! carbon dioxide (mmr)
     |    n2   (levd0,lond0,latd0,itc),      ! updated n2 (mmr) from compn2
     |    so   (levd0,lond0,latd0,itp),      ! previous SO (at n) on history
     |    so_nm (levd0,lond0,latd0,itp),     ! previous SO (at n-1) on history
     |    so   (levd0,lond0,latd0,itc),      ! new SO (at n) on history
     |    so_nm (levd0,lond0,latd0,itc),     ! new SO (at n-1) on history
     |    1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('driver after minor_so')")
! --------------------------------------------------------------------------
! Advance so2 (minor): 14-inputs to call
!     subroutine comp_so2(tn,o1,co,co2,n2,o2,so,so2,barm,
!    |  lev0,lev1,lon0,lon1,lat)
! --------------------------------------------------------------------------
       do lat=lat0,lat1
        if (debug) write(6,"('driver calling comp_so2: lat=',i3)") lat
        call comp_so2(
     |    tn   (levd0,lond0,lat,itc), ! neutral temperature
     |    o1   (levd0,lond0,lat,itc), ! atomic oxygen (mmr)
     |    co   (levd0,lond0,lat,itc), ! carbon momoxide (mmr)
     |    co2  (levd0,lond0,lat),     ! carbon dioxide (mmr)
     |    n2   (levd0,lond0,lat,itc), ! updated psi n2 from compn2
     |    o2   (levd0,lond0,lat,itc), ! updated psi o2 from compo2
!    |    so   (levd0,lond0,lat,itp), ! so (old, from history since yet call comp_so)
     |    so   (levd0,lond0,lat,itc), ! so (new, on history since call comp_so)
     |    so2  (levd0,lond0,lat,itp), ! so2 (old, from history)
     |    barm (levd0,lond0,lat),     ! mean molecular weight
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('driver after comp_so2')")
! --------------------------------------------------------------------------
!     subroutine minor_so2(tn,difk,o1,co,co2,n2,so2,so2_nm1,so2_out,so2_nm1_out,
!    |  lev0,lev1,lon0,lon1,lat0,lat1)
! --------------------------------------------------------------------------
        call minor_so2(
     |    tn   (levd0,lond0,latd0,itc),      ! neutral temperature
     |    difk (levd0,lond0,latd0),          ! eddy diffusion coefficient (eddyflds)
     |    o1   (levd0,lond0,latd0,itc),      ! atomic oxygen (mmr)
     |    co   (levd0,lond0,latd0,itc),      ! carbon momoxide (mmr)
     |    co2  (levd0,lond0,latd0),          ! carbon dioxide (mmr)
     |    n2   (levd0,lond0,latd0,itc),      ! updated n2 (mmr) from compn2
     |    so2  (levd0,lond0,latd0,itp),      ! previous SO2 (at n) on history
     |    so2_nm (levd0,lond0,latd0,itp),    ! previous SO2 (at n-1) on history
     |    so2  (levd0,lond0,latd0,itc),      ! new SO2 (at n) on history
     |    so2_nm (levd0,lond0,latd0,itc),    ! new SO2 (at n-1) on history
     |    1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('driver after minor_so2')")
!
! -----------------------------------------------------------------
!
! so2 and so are redundant in time until the comp and minor subs are ready.
! (so2 and so were initialized at itp and istep=1 in chapman.F)
!
!     so2(:,:,:,itc) = so2(:,:,:,itp)
!     so2_nm(:,:,:,itc) = so2_nm(:,:,:,itp)
!     so (:,:,:,itc) = so (:,:,:,itp)
!     so_nm (:,:,:,itc) = so_nm (:,:,:,itp)
!
! Sources and sinks for major species composition OH:
!     subroutine comp_oh(tn,barm,o1,co,co2,n2,o2,op,co2p,o2p,ne,
!    |  lev0,lev1,lon0,lon1,lat)
! -----------------------------------------------------------------
      do lat=lat0,lat1
        if (debug) write(6,"('driver calling comp_oh: lat=',i3)") lat
        call comp_oh(
     |    tn   (levd0,lond0,lat,itp), ! neutral temperature
     |    barm (levd0,lond0,lat),     ! mean molecular weight
     |    o1   (levd0,lond0,lat,itp), ! atomic oxygen (mmr)
     |    co   (levd0,lond0,lat,itp), ! carbon momoxide (mmr)
     |    co2  (levd0,lond0,lat),     ! carbon dioxide (mmr)
     |    n2   (levd0,lond0,lat,itc), ! updated psi n2 from compn2
     |    o2   (levd0,lond0,lat,itp), ! molecular oxygen (mmr)
     |    op   (levd0,lond0,lat,itc), ! O+ (on history)
     |    co2p (levd0,lond0,lat),     ! CO2+ (diagnostic)
     |    o2p  (levd0,lond0,lat,itc), ! O2+ (on history)
     |    ne   (levd0,lond0,lat,itc), ! electron density (on history)
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('driver after comp_oh')")
! --------------------------------------------------------------------------
!     do lat=lat0,lat1
!       call addfld('Z_end',' ',' ',z(:,lon0:lon1,lat,itc),
!    |    'lev',1,nlevp1,'lon',lon0,lon1,lat)
!     enddo
      end subroutine driver
