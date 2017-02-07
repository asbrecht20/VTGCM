!
      subroutine duv(tn,tn_upd,tn_nm,un,vn,un_nm,vn_nm,w_upd,barm,
     |  xmue,z,hdu,hdv,km,un_upd,unm_upd,vn_upd,vnm_upd,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
! -----------------------------------------------------------------------------
! Advance neutral velocities for current time step:
! 9/21/06 btf: lamdas and ion drifts are made local rather than in
!              arg list, and are zeroed at near beginning of this routine.
! 11/14/06 swb: add RF and RSZ prescribed flow code from vtgcm2
!               Temporary formulation until gwdd.F can be incorporated
!               Symmetric tdrag (start); gwdd.F provides RSZ flow (later)
! 12/07/06 swb: add xmue to arg list (from eddyflds)
! 12/22/06 swb: tscale = 0.8E-04 (1/sec) weakened for smin conditions
! 1/16/07  bf:  mods for new_lbc_uv (ubound, vbound)
! 01/23/07 swb: add RF and RSZ prescribed flow code from vtgcm2
!               Temporary formulation until gwdd.F can be incorporated
!               Asymmetric drag (upgrade) for u-momentum (only). U(RSZ).
! 09/12/07  swb: tscale = 1.0E-04 (1/sec) weakened for smod conditions
! 01/30/08  swb: tscale = 1.2E-04 (1/sec), zpbk = +2.0  (best)
! 02/13/08  swb: tscale = 1.3E-04 (1/sec), zpbk = +2.0  (best)
! 02/14/08  swb: tscale = 1.2E-04 (1/sec), zpbk = +2.0  (best)
! 03/07/08  asb: tscale = 0.1E-05 (1/sec), zpbk = +2.0  (extreme)
! 10/15/08  asb: tscale = 0.7E-04 (1/sec), zpbk = +2.0  (overall best)
! 02/20/09  swb: URSZ (Explicit drag term) added to RHS of U-mom Equation
!              :  RHS = RHS - tdrag(k)*ursz*cs(lat)*vrt(k)
!              :  ursz = +7500.0 cm/sec  (at zpbk = 2 and above, VRT = profile)
!              :  ursz = +7500.0 cm/sec  (at zpbk = 1 and above, VRT = profile)
! 02/23/09  swb: URSZ (Explicit drag term) added to RHS of U-mom Equation
!              :  ursz = +7500.0 cm/sec  (at zpbk = 2 and above, VRT2 = 1.0)
!              :  ursz = +7500.0 cm/sec (at zpbk = 2 and above, VRT2 = SQRT(VRT)
!              :  ursz = +7500.0 cm/sec  (at zrszf = 0, VRT2 = SQRT(1-expo)
! 02/24/09  swb: URSZ (Explicit drag term) added to RHS of U-mom Equation
!              :  ursz = +7500.0 cm/sec  (at zrszf = -1.5, VRT2 = SQRT(1-expo)
! 05/14/14  swb:  temporary U=V=0.0 set for diagnostic (1-D) profile testing
! 05/28/14  swb:  interative winds back on
! 06/02/14  swb:  temporary U=V=0.0 set for diagnostic (1-D) profile testing
! 06/02/14  swb:  all self-consistent winds back on;  weak tdrag (fast winds)
! 06/13/14  swb:  all self-consistent winds back on;  tscale = 1.2e-4 (nominal winds)
! 06/17/14  swb:  all self-consistent winds back on;  tscale = 1.25e-5 (fast winds)
! 06/17/14  swb:  all self-consistent winds back on;  tscale = 1.0e-4 (nominal winds)
! -----------------------------------------------------------------------------
!
      use params_module,only: nlonp4,dz,nlat,spval,zibot,nlevp1
      use init_module,only: iter
      use lbc,only: u_lbc,v_lbc
      use cons_module,only: freq_semidi,dt,dtx2inv,expz,p0,gask,grav,
     |  cor,tanphi=>tn,re,kut_5,kut,dtsmooth_div2,dtsmooth,cs
      use fields_module,only: tlbc,ulbc,vlbc,ulbc_nm,vlbc_nm
      use input_module,only: oxvgcm_ncfile,oxvgcm_zm
      use oxvgcm,only: oxvgcm_nlev,oxvgcm_nudge,u_oxvgcm,v_oxvgcm,
     |  uzm_oxvgcm,vzm_oxvgcm
      use addfld_module,only: addfld
#ifdef MPI
      use mpi_module,only: mp_periodic_f3d
#endif
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
!
! Inputs at full subdomain:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn,      ! neutral temperature (deg K)
     |  tn_upd,  ! updated neutral temperature (from dt.F)
     |  tn_nm,   ! neutral temperature from time n-1
     |  un_nm,   ! neutral zonal velocity at time n-1
     |  vn_nm,   ! neutral meridional velocity at time n-1
     |  w_upd,   ! updated vertical velocity [/sec](swdot.F)
     |  barm,    ! mean molecular weight
     |  xmue,    ! eddy viscosity coefficient (1/sec) (eddyflds)
     |  z,       ! geopotential height
     |  hdu,     ! horizontal diffusion of U (hdif3 in hdif.F)
     |  hdv,     ! horizontal diffusion of V (hdif3 in hdif.F)
     |  km       ! molecular viscosity (cpktkm.F)
!
! un,vn are intent(inout) so they can be nudged at beginning:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(inout)::
     |  un,      ! neutral zonal velocity
     |  vn       ! neutral meridional velocity
!
! Outputs at full subdomain:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out)::
     |  un_upd,  ! updated zonal velocity
     |  unm_upd, ! updated zonal velocity at time n-1
     |  vn_upd,  ! updated meridional velocity
     |  vnm_upd  ! updated meridional velocity at time n-1
!
! Local:
      integer :: k,i,j,lat,lonbeg,lonend,kutt(nlat)
      real :: 
     |  unlbc(lon0:lon1), ! un lower boundary (t1)
     |  vnlbc(lon0:lon1)  ! vn lower boundary (t2)
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  g,           ! g**2*(kM+XMUE)*MBAR*/(P0*R*T*DS**2)
     |  dwdz,        ! exp(-s)*w(k+1/2)/(2.*Ds)
     |  tni,         ! tn at interfaces
     |  hsc,         ! scale height [cm]
     |  lxx,lyy,lxy,lyx,ui,vi
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2) ::
     |  w_updm,      ! updated vertical velocity [m/sec] (swdot.F)
     |  dragtermdcs, ! -tdrag*(un_upd-2*vrt2*cs)/100 [m/sec2]
     |  dragtermcs     ! -tdrag*(un_upd-vrt2*cs)/100 [m/sec2]
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  advec_un,    ! zonal advection (output of sub advec)
     |  advec_vn,    ! meridional advection (output of sub advec)
     |  zl,zp,       ! output from glp (and dldp)
     |  unm_smooth,  ! un at time n-1, smoothed
     |  vnm_smooth,  ! vn at time n-1, smoothed
     |  rtxmue       ! sqrt(xmue)
!
! Zonal means for oxvgcm nudging:
      real,dimension(lev0:lev1,lat0:lat1) :: uzm,vzm ! zonal means for nudging
      real,dimension(lon0:lon1) :: uzm_pert,vzm_pert ! perturations about mean
!
!     real :: rtxmue(lev0:lev1) ! sqrt(xmue)
!!    real :: tdrag(lev0:lev1), vrt(lev0:lev1), vrt2(lev0:lev1)  ! RFFormulation
      real :: tdrag(lev0:lev1), vrt(lev0:lev1)   ! RFFormulation
      real :: vrt2(nlevp1)  ! RFFormulation

!     real :: ursz(lev0:lev1)
      real :: tscale, zpbk, szp, zrszf,ursz  !  R.FRICTION Formulation
!
! For diagnostic, debug:
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  ss_un,ss_vn,  ! also for diagnostics
     |  pp11,pp12,pp21,pp22, qq11,qq12,qq21,qq22, rr11,rr12,rr21,rr22
      real,dimension(lon0:lon1,lat0:lat1) ::
     |  unlbc_ij,     ! un lbc , for diagnostic ij
     |  vnlbc_ij      ! vn lbc , for diagnostic ij
!
! diagnostics for zonal mean nudging:
      real :: uzm_vtgcm_diag(lev0:lev1,lon0:lon1,lat0:lat1)
      real :: uzm_oxvgcm_diag(lev0:lev1,lon0:lon1,lat0:lat1)
      real :: vzm_vtgcm_diag(lev0:lev1,lon0:lon1,lat0:lat1)
      real :: vzm_oxvgcm_diag(lev0:lev1,lon0:lon1,lat0:lat1)
!
! pp, qq, rr, beta AND gamma are (2 x 2) matrices defined for 
!   k = 3/2, K-1/2, 1 and i = 3, IMAX+2, 1
      real,dimension(2,2,lev0:lev1,lon0:lon1) :: 
     |  pp,qq,rr,beta,gamma 
!
! ss, xx and yy are two component vectors similarly defined.
      real,dimension(2,lev0:lev1,lon0:lon1) :: 
     |  ss,xx,yy
!
#include<defs.h>
!
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = nlonp4-2
!
! Calculate zonal means for oxvgcm nudging:
      if (len_trim(oxvgcm_ncfile) > 0.and.oxvgcm_zm > 0) then
        call calc_uvzm(un(:,lon0:lon1,lat0:lat1),uzm,
     |                 vn(:,lon0:lon1,lat0:lat1),vzm,
     |    lev0,lev1,lon0,lon1,lat0,lat1)
!
! Plot zonal means (diags are redundant in longitude):
        do j=lat0,lat1
          do i=lon0,lon1
            uzm_vtgcm_diag(:,i,j) = uzm(:,j)          ! redundant in lon
            uzm_oxvgcm_diag(:,i,j) = uzm_oxvgcm(j,:) ! redundant in lon
            vzm_vtgcm_diag(:,i,j) = vzm(:,j)          ! redundant in lon
            vzm_oxvgcm_diag(:,i,j) = vzm_oxvgcm(j,:) ! redundant in lon
          enddo
          call addfld('uzm_vtgcm',' ',' ',uzm_vtgcm_diag(:,:,j),
     |      'lev',lev0,lev1,'lon',lon0,lon1,j)
          call addfld('uzm_oxvgcm',' ',' ',uzm_oxvgcm_diag(:,:,j),
     |      'lev',lev0,lev1,'lon',lon0,lon1,j)
          call addfld('vzm_vtgcm',' ',' ',vzm_vtgcm_diag(:,:,j),
     |      'lev',lev0,lev1,'lon',lon0,lon1,j)
          call addfld('vzm_oxvgcm',' ',' ',vzm_oxvgcm_diag(:,:,j),
     |      'lev',lev0,lev1,'lon',lon0,lon1,j)
        enddo
      endif
!
! R. Friction constants: vtgcm2 from 1997 Benchmark simulations
! zpbk, ursz (cm/sec), tscale (sec-1)
! ------------------------
!     ursz = 7500.0
!     ursz = 9000.0
!     ursz = 0.0
      ursz = 1.0
!     zrszf = 0
!     zrszf = -0.5
!     zrszf = -1  !Best Case for Aussois
!     zrszf = 1.
!     zrszf = 3.0
!     zrszf = 4.0
!     zrszf = -1.5
!     zrszf = -2.5
!     zrszf = -3.
! ------------------------
!     zpbk = 2.0
!     zpbk = 1.5  !Best Case for Aussois
      zpbk = -1.5
!     zpbk = -2.5
!     zpbk = 0
!     zpbk = -1
!     zpbk = -3
!     zpbk = -4
!     zpbk = 0.5
! ------------------------
!     Extreme case (230K on nightside @ 100km)
!     tscale = 0.05E-04
!     tscale = 0.01E-04
!     Normal case
!     tscale = 1.2E-04
!     tscale = 1.08E-04  !used to experiment with NO*
!     tscale = 1.3E-04
!     tscale = 0.9E-04
!     tscale = 3.0E-04   !Looking for Min wind case Night Temp ~105km = 175K
!     tscale = 0.1E-04   !Max wind conditions Night Temp ~100km = 180K
!     tscale = 0.7E-04   !Norm Case Night Temp @100km = 200K
!     tscale = 0.5E-04   ! Best case for O2 nightglow position
!     tscale = 0.
!     tscale = 1.25E-05  ! Min case for O2 nightglow position (fast winds)
!     tscale = 1.2E-04   ! VTGCM old case  (best wave drag for Tnite) 
!     tscale = 2.0E-04   ! Max case for O2 nightglow position (slow winds)
!     tscale = 1.4E-04   ! Minimal case Night Temp @100 = 190K
      tscale = 1.0E-04   ! Best overall case
!     tscale = 1.8E-04
!     tscale = 1.5E-04
!     tscale = 2.0E-04   !Min wind conditions  Night Tem @100 = 180K
! ------------------------
!
! Zero out lamdas and ion drift velocities:
      lxx = 0.
      lyy = 0.
      lxy = 0.
      lyx = 0.
      ui = 0.
      vi = 0.
! Implicit Portion
! Setup Rayleigh Friction formulation (versus altitude)
! Taken from the vtgcm2 duv.F code (1997)
! tdrag ----------------------------------
      do k=lev0,lev1-1
	    szp = zibot + .5*dz + (k-1)*dz
	 if(szp >= zpbk) then
	    vrt(k) = 1.0
         else
!           vrt(k) = 2.*exp((szp-zpbk)/2.)/(1.+exp((szp-zpbk)/2.))
!           vrt(k) = 2.*sqrt(exp((szp-zpbk)/2.)/(2.+exp((szp-zpbk)/2.)))
 	    vrt(k) = exp((szp-zpbk)/2.)
         endif
	    tdrag(k) = vrt(k)*tscale
      enddo ! k=lev0,lev1
! ursz profile ----------------------------------
!
!!    do k=lev0,lev1-1
!!          szp = zibot + .5*dz + (k-1)*dz
!!       if(szp >= zrszf) then
!!          vrt2(k) = 1.0
!!       else
!           vrt2(k) = sqrt(exp((szp-zrszf)/2.))
!!          vrt2(k) = exp((szp-zrszf)/2.)
!           vrt2(k) = exp((szp-zrszf))
!****vrt2*ursz prescribed ****
#if (NLEV==68)
      vrt2 = (/
     |  1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
     |  1., 1., 1., 1., 1., 1., 1., 1.,
     |  1., 1., 1., 1., 1., 100.,
     |  2226.99, 2849.53, 3646.08, 4665.30, 5969.44, 7500.,
     |  7500., 7500., 7500., 7500., 7500., 7500., 7500.,
     |  7500., 7500., 7500., 7500., 7500., 7500.,
     |  7500., 7500., 7500., 7500., 7500., 7500., 7500.,
     |  7500., 7500., 7500., 7500., 7500., 7500.,
     |  7500., 7500., 7500., 7500., 7500., 7500., 7500.,
     |  7500., 7500., 7500., 7500., 7500./)
#elif (NLEV==136)
      vrt2 = (/
     |  1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 
     |  1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 
     |  1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
     |  1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
     |  1., 1., 1., 1., 100., 1163.50, 2226.99, 
     |  2538.26, 2849.53, 3247.81, 3646.08, 4155.7,
     |  4665.30, 5317.37, 5969.44, 6734.72, 7500.,
     |  7500.,7500., 7500., 7500., 7500., 7500., 7500., 
     |  7500., 7500., 7500., 7500., 7500., 7500.,
     |  7500., 7500., 7500., 7500., 7500., 7500., 7500., 
     |  7500., 7500., 7500., 7500., 7500., 7500., 
     |  7500., 7500., 7500., 7500., 7500., 7500., 7500., 
     |  7500., 7500., 7500., 7500., 7500.,
     |  7500., 7500., 7500., 7500., 7500., 7500., 7500., 
     |  7500., 7500., 7500., 7500., 7500., 7500., 7500.,
     |  7500., 7500., 7500., 7500., 7500., 7500., 7500., 
     |  7500., 7500., 7500., 7500., 7500., 7500., 
     |  7500., 7500., 7500., 7500., 7500., 7500., 7500., 
     |  7500., 7500., 7500., 7500./)
#else
 write(6,"('ERROR: unsupported NLEV = ',i6)") nlev
#endif
!!       endif
!!    enddo ! k=lev0,lev1
! ----------------------------------------------
!
! First latitude scan:
      do lat=lat0,lat1
!   
! Nudge with Oxford VGCM (u,v_oxvgcm(nlonp4,nlat,oxvgcm_nlev)). See lbc.F.
! oxvgcm_nudge==1.0 at k==1, decreasing upward to 0 at k==oxvgcm_nlev.
! When nudge is 1, then we use 100% oxvgcm, and when nudge is 0, then
!   we use 100% vtgcm. Bottom boundary is already 100% oxvgcm, so start
!   at k==2.
! Nudge with zonal means of oxvgcm if namelist parameter oxvgcm_zm > 0
! (zonal means of u,v were calculated in call to calc_uvzm above)
!    
        if (len_trim(oxvgcm_ncfile) > 0) then

          call addfld('u_prenudge',' ',' ',un(:,lon0:lon1,lat),
     |      'lev',lev0,lev1,'lon',lon0,lon1,lat)
          call addfld('v_prenudge',' ',' ',vn(:,lon0:lon1,lat),
     |      'lev',lev0,lev1,'lon',lon0,lon1,lat)

          if (oxvgcm_zm > 0) then   ! nudge zonal mean
            do k=2,oxvgcm_nlev
!
! Save perturbations about the mean:
              do i=lon0,lon1
                uzm_pert(i) = un(k,i,lat)-uzm(k,lat)
                vzm_pert(i) = vn(k,i,lat)-vzm(k,lat)
              enddo
!
! Nudge zonal mean uzm,vzm (from sub calc_uvzm) with oxvgcm zonal means:
              uzm(k,lat) = uzm(k,lat)*(1.-oxvgcm_nudge(k))+
     |          uzm_oxvgcm(lat,k)*oxvgcm_nudge(k)
              vzm(k,lat) = vzm(k,lat)*(1.-oxvgcm_nudge(k))+
     |          vzm_oxvgcm(lat,k)*oxvgcm_nudge(k)
!
! Update u,v w/ nudged zm and restore perturbations about the mean:
              do i=lon0,lon1
                un(k,i,lat) = uzm(k,lat)+uzm_pert(i)
                vn(k,i,lat) = vzm(k,lat)+vzm_pert(i)
              enddo
            enddo ! k=2,oxvgcm_nlev

            call addfld('u_zmnudge',' ',' ',un(:,lon0:lon1,lat),
     |        'lev',lev0,lev1,'lon',lon0,lon1,lat)
            call addfld('v_zmnudge',' ',' ',vn(:,lon0:lon1,lat),
     |        'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Normal nudging with u,v_oxvgcm:
          else
            do k=2,oxvgcm_nlev
              do i=lon0,lon1
                un(k,i,lat) = un(k,i,lat)*(1.-oxvgcm_nudge(k))+
     |            u_oxvgcm(i,lat,k)*oxvgcm_nudge(k)
                vn(k,i,lat) = vn(k,i,lat)*(1.-oxvgcm_nudge(k))+
     |            v_oxvgcm(i,lat,k)*oxvgcm_nudge(k)
              enddo
            enddo
!           write(6,"('duv nudging un,vn with oxvgcm: lat=',i3,
!    |        ' oxvgcm_nudge=',/,(6e12.4))") lat,oxvgcm_nudge
          endif ! nudging w/ zonal means
        endif ! oxvgcm nudging
!
! Horizontal advection:
!     subroutine advec(f,hadvec,lev0,lev1,lon0,lon1,lat)
!     real,dimension(lev0:lev1,lon0-2:lon1+2,lat-2:lat+2),intent(in) :: 
!    |  f ! input field with ghost cells for finite differencing
!     real,dimension(lev0:lev1,lon0:lon1),intent(out) :: hadvec
!      
        call advec(un(:,:,lat-2:lat+2),advec_un(:,:,lat), 
     |    lev0,lev1,lon0,lon1,lat)
        call advec(vn(:,:,lat-2:lat+2),advec_vn(:,:,lat), 
     |    lev0,lev1,lon0,lon1,lat)
!
!       call addfld('ADVEC_U0',' ',' ',advec_un(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('ADVEC_V0',' ',' ',advec_vn(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! End first latitude scan:
      enddo ! lat=lat0,lat1
!
! Horizontal pressure forcing output of glp is in zl,zp. Pass full 
! subdomains for derivatives in lat and lon (see sub dldp, called by glp):
!
!     do lat=lat0,lat1
!       call addfld('DUV_Z',' ',' ',z(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     enddo

      call glp(tn,tn_upd,tn_nm,barm,z,zl,zp,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Second latitude scan:
      do lat=lat0,lat1

!       call addfld('ZL',' ',' ',zl(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('ZP',' ',' ',zp(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)

!       do k = lev0,lev1 - 1
!         write(6,"('duv: lat=',i3,' k=',i3,' zl(k,:,lat)=',/,
!    |      (6e12.4))") lat,k,zl(k,:,lat)
!         write(6,"('  zp(k,:,lat)=',/,(6e12.4))") zp(k,:,lat)
!       enddo

        do i=lon0,lon1
          do k=lev0,lev1-1
            advec_un(k,i,lat) = advec_un(k,i,lat)+
     |        zl(k,i,lat)-hdu(k,i,lat)
            advec_vn(k,i,lat) = advec_vn(k,i,lat)+
     |        zp(k,i,lat)-hdv(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       call addfld('ADVEC_U1',' ',' ',advec_un(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('ADVEC_V1',' ',' ',advec_vn(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! End second lat scan:
      enddo ! lat=lat0,lat1
!
! Smooth un,vn at time n-1 with 2-part shapiro smoother.
! Input un_nm,vn_nm at full subdomain.
! Output in unm_smooth,vnm_smooth at subdomain.
!
      call smooth(un_nm,unm_smooth,lev0,lev1,lon0,lon1,lat0,lat1,0)
      call smooth(vn_nm,vnm_smooth,lev0,lev1,lon0,lon1,lat0,lat1,0)
!
! Start third latitude scan:
      do lat=lat0,lat1

!       call addfld('UNSMOOTH',' ',' ',unm_smooth(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('VNSMOOTH',' ',' ',vnm_smooth(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! ------------------------------------------------------------------
! U,V bottom boundary:
! ------------------------------------------------------------------
!
!
! Lower boundaries u_lbc,v_lbc were calculated by sub tuvz_lbc (lbc.F)
        unlbc(:) = u_lbc(lon0:lon1,lat)
        vnlbc(:) = v_lbc(lon0:lon1,lat)
!
! Define unlbc_ij for diagnostic:
        do i=lon0,lon1
          unlbc_ij(i,lat) = unlbc(i)
          vnlbc_ij(i,lat) = vnlbc(i)
        enddo ! i=lon0,lon1
!
! Save 2d lon,lat boundaries to secondary history:
        call addfld('UN_LBC','UN Lower boundary','cm/s',unlbc_ij,
     |    'lon',lon0,lon1,'lat',lat0,lat1,0)
        call addfld('VN_LBC','VN Lower boundary','cm/s',vnlbc_ij,
     |    'lon',lon0,lon1,'lat',lat0,lat1,0)
!
! ------------------------------------------------------------------
!  Add Explicit Drag Term Here:  RHS term (ss(1,k,i))
!  .....tdrag(k)*ursz*cs(lat)*vrt(k)
!  .....tdrag(k)*ursz*cs(lat)*1.0
!
        do i=lonbeg,lonend
          do k=lev0,lev1-1
! -----------------------------------------------------------------
!**No Asymmetrical Drag**
!           advec_un(k,i,lat) = advec_un(k,i,lat)-
!    |        dtx2inv*unm_smooth(k,i,lat)
!           advec_un(k,i,lat) = advec_un(k,i,lat)-
!    |        dtx2inv*unm_smooth(k,i,lat)-
!    |        tdrag(k)*ursz*cs(lat)*vrt(k)
!         advec_un(k,i,lat) = advec_un(k,i,lat)-
!    |        dtx2inv*unm_smooth(k,i,lat)-
!    |        tdrag(k)*ursz*cs(lat)*1.0
!**Asymmetrical Drag**
          advec_un(k,i,lat) = advec_un(k,i,lat)-
     |        dtx2inv*unm_smooth(k,i,lat)-
     |        tdrag(k)*ursz*cs(lat)*vrt2(k)
!         advec_un(k,i,lat) = advec_un(k,i,lat)-
!    |        dtx2inv*unm_smooth(k,i,lat)-
!    |        tdrag(k)*ursz*cs(lat)
! -----------------------------------------------------------------
            advec_vn(k,i,lat) = advec_vn(k,i,lat)-
     |        dtx2inv*vnm_smooth(k,i,lat)
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1

!       call addfld('ADVEC_U2',' ',' ',advec_un(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('ADVEC_V2',' ',' ',advec_vn(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)

        do i=lon0,lon1
          do k=lev0,lev1-1
            ss(1,k,i) = expz(k)*(.25*((lxx(k,i)+lxx(k+1,i))*
     |        (ui(k,i)+ui(k+1,i))+(lxy(k,i)+lxy(k+1,i))*
     |        (vi(k,i)+vi(k+1,i)))-advec_un(k,i,lat))
            ss(2,k,i) = expz(k)*(.25*((lyy(k,i)+lyy(k+1,i))*
     |        (vi(k,i)+vi(k+1,i))-(lyx(k,i)+lyx(k+1,i))*
     |        (ui(k,i)+ui(k+1,i)))-advec_vn(k,i,lat))
            ss_un(k,i) = ss(1,k,i)
            ss_vn(k,i) = ss(2,k,i)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       call addfld('SS1_UN',' ',' ',ss_un,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('SS1_VN',' ',' ',ss_vn,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)

!
! xmue = eddy viscosity (eddyflds.F). rtxmue = sqrt(xmue)
        do i=lon0,lon1
          do k=lev0+1,lev1-1          
           rtxmue(k,i,lat) = sqrt(xmue(k-1,i,lat)*xmue(k,i,lat))
          enddo ! k=lev0+1,lev1-1          
!
! Boundaries:
          rtxmue(lev0,i,lat) = sqrt(xmue(lev0,i,lat)**3/
     |                         xmue(lev0+1,i,lat))
          rtxmue(lev1,i,lat) = sqrt(xmue(lev1-1,i,lat)**3/
     |                         xmue(lev1-2,i,lat))
        enddo ! i=lon0,lon1
!
! tni = tn at interfaces:
        do i=lon0,lon1
          do k=lev0+1,lev1-1          
            tni(k,i) = .5*(tn(k-1,i,lat)+tn(k,i,lat))
          enddo ! k=lev0+1,lev1-1          
          tni(lev0,i) = tlbc(i,lat)      ! bottom boundary in tlbc
          tni(lev1,i) = tn(lev1-1,i,lat) ! lev1 <- lev1-1
!
! G = g*KM/(P0*H*Ds**2) = g**2*(kM+XMUE)*MBAR*/(P0*R*T*DS**2)
          do k=lev0,lev1
            g(k,i) = grav**2*(km(k,i,lat)+rtxmue(k,i,lat))*
     |        barm(k,i,lat)/(p0*gask*tni(k,i)*dz**2)
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1

!       call addfld('TNI'  ,' ',' ',tni,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('DUV_G',' ',' ',g  ,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! dwdz = exp(-s)*w(k+1/2)/(2.*Ds)
        do i=lon0,lon1
          do k=lev0,lev1-1
            dwdz(k,i) = 0.5*expz(k)*(w_upd(k,i,lat)+w_upd(k+1,i,lat))/
     |        (2.*dz)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       call addfld('DWDZ',' ',' ',dwdz,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! P = -(G(k-1/2,n) + exp(-s)*w(k,n)/   |1., 0.|
!                            (2.*Ds)) *|0., 1.|
!                                      |_    _|
!
! R = -(G(k+1/2,n) - exp(-s)*w(k,n)/   |1., 0.|
!                            (2.*Ds)) *|0., 1.|
!                                      |_    _|
!  for k = 3/2, K-1/2, 1
!
        do i=lon0,lon1
          do k=lev0,lev1-1
            pp(1,1,k,i) = -(g(k  ,i)+dwdz(k,i))
            rr(1,1,k,i) = -(g(k+1,i)-dwdz(k,i))
            pp(2,2,k,i) = pp(1,1,k,i)
            rr(2,2,k,i) = rr(1,1,k,i)
            pp(1,2,k,i) = 0.
            rr(1,2,k,i) = 0.
            pp(2,1,k,i) = 0.
            rr(2,1,k,i) = 0.
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
!
! Now calculate Q from (5) for k = 3/2, K-1/2, 1
!  Q(1,1,k) = G(k-1/2) + G(k+1/2) + exp(-s)*(1./(2.*Dt) + lamxx(k))
!    and
!  Q(2,2,k) = G(k-1/2) + G(k+1/2) + exp(-s)*(1./(2.*Dt) + lamyy(k))
!
!     Location of RF formulation for tdrag (from vtgcm2)
!     (11/14/06) SWB: Patterned after jtgcm16 coding
!                   : Symmetric tdrag (only) to start
!
        do i=lon0,lon1
          do k=lev0,lev1-1
            qq(1,1,k,i) = g(k,i) + g(k+1,i)
            qq(2,2,k,i) = qq(1,1,k,i)
! -------------------------------------------------------------------
! Asymmetric (only), with LAT and ALT variation :  u-momentum  (NG)
!           qq(1,1,k,i) = qq(1,1,k,i) +
!    |        expz(k)*(dtx2inv+.5*(lxx(k,i)+lxx(k+1,i))+
!    |        tdrag(k)*(1.-0.15*cs(lat)*vrt(k)))
! -------------------------------------------------------------------
! Symmetric (always) :  u-momentum
            qq(1,1,k,i) = qq(1,1,k,i) +
!    |        expz(k)*(dtx2inv+.5*(lxx(k,i)+lxx(k+1,i)))
     |        expz(k)*(dtx2inv+.5*(lxx(k,i)+lxx(k+1,i))+tdrag(k))
! Symmetric (always): v-momentum
            qq(2,2,k,i) = qq(2,2,k,i) +
!    |        expz(k)*(dtx2inv+.5*(lyy(k,i)+lyy(k+1,i)))
     |        expz(k)*(dtx2inv+.5*(lyy(k,i)+lyy(k+1,i))+tdrag(k))
! -------------------------------------------------------------------
!
!     Q(1,2) = -exp(-s) * (f+u/r*tan(theta)-lamxy)
!      and 
!     Q(2,1) = exp(-s) * (f+u/r*tan(theta)-lamyx)
!
            qq(1,2,k,i) = cor(lat)+un(k,i,lat)/re*tanphi(lat)
            qq(2,1,k,i) = qq(1,2,k,i)
            qq(1,2,k,i) = -expz(k)*(qq(1,2,k,i)-
     |        .5*(lxy(k,i)+lxy(k+1,i)))
            qq(2,1,k,i) =  expz(k)*(qq(2,1,k,i)-
     |        .5*(lyx(k,i)+lyx(k+1,i)))
          enddo ! k=lev0,lev1-1
!
! Boundaries at k=3/2 and k=k-1/2
! Lower boundary:
!     Q(3/2) = Q(3/2) - P(3/2)
!     S(3/2) = S(3/2) -2.*P(3/2)*Vb
!     P(3/2) = 0.
!
          qq(1,1,lev0,i) = qq(1,1,lev0,i)-pp(1,1,lev0,i)
          qq(1,2,lev0,i) = qq(1,2,lev0,i)-pp(1,2,lev0,i)
          qq(2,1,lev0,i) = qq(2,1,lev0,i)-pp(1,2,lev0,i)
          qq(2,2,lev0,i) = qq(2,2,lev0,i)-pp(2,2,lev0,i)
          ss(1,lev0,i) = ss(1,lev0,i)-2.*(pp(1,1,lev0,i)*unlbc(i)+
     |                                    pp(1,2,lev0,i)*vnlbc(i))
          ss(2,lev0,i) = ss(2,lev0,i)-2.*(pp(2,1,lev0,i)*unlbc(i)+
     |                                    pp(2,2,lev0,i)*vnlbc(i))
          ss_un(lev0,i) = ss(1,lev0,i)
          ss_vn(lev0,i) = ss(2,lev0,i)
          pp(1,1,lev0,i) = 0.
          pp(1,2,lev0,i) = 0.
          pp(2,1,lev0,i) = 0.
          pp(2,2,lev0,i) = 0.
!
! Upper boundary:
!     Q(K-1/2) = Q(K-1/2) + R(K-1/2)
!     R(K-1/2) = 0.
!
          qq(1,1,lev1-1,i) = qq(1,1,lev1-1,i)+rr(1,1,lev1-1,i)
          qq(1,2,lev1-1,i) = qq(1,2,lev1-1,i)+rr(1,2,lev1-1,i)
          qq(2,1,lev1-1,i) = qq(2,1,lev1-1,i)+rr(2,1,lev1-1,i)
          qq(2,2,lev1-1,i) = qq(2,2,lev1-1,i)+rr(2,2,lev1-1,i)
          rr(1,1,lev1-1,i) = 0.
          rr(1,2,lev1-1,i) = 0.
          rr(2,1,lev1-1,i) = 0.
          rr(2,2,lev1-1,i) = 0.
        enddo ! i=lon0,lon1
!
! diags:
!       pp11(:,:) = pp(1,1,:,:) ; pp12(:,:) = pp(1,2,:,:)
!       pp21(:,:) = pp(2,1,:,:) ; pp22(:,:) = pp(2,2,:,:)
!       qq11(:,:) = qq(1,1,:,:) ; qq12(:,:) = qq(1,2,:,:)
!       qq21(:,:) = qq(2,1,:,:) ; qq22(:,:) = qq(2,2,:,:)
!       rr11(:,:) = rr(1,1,:,:) ; rr12(:,:) = rr(1,2,:,:)
!       rr21(:,:) = rr(2,1,:,:) ; rr22(:,:) = rr(2,2,:,:)
!       call addfld('PP11',' ',' ',pp11,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('PP12',' ',' ',pp12,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('PP21',' ',' ',pp21,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('PP22',' ',' ',pp22,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('QQ11',' ',' ',qq11,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('QQ12',' ',' ',qq12,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('QQ21',' ',' ',qq21,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('QQ22',' ',' ',qq22,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('RR11',' ',' ',rr11,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('RR12',' ',' ',rr12,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('RR21',' ',' ',rr21,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('RR22',' ',' ',rr22,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('SS2_UN',' ',' ',ss_un,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('SS2_VN',' ',' ',ss_vn,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Invoke tridiagonal (2x2) block matrix solver for U and V: 
! pp,qq,rr,ss are input, beta,gamma,yy,xx are output.
! (final solution is in xx)
!
        call blktri(pp,qq,rr,ss, beta,gamma,yy,xx,
     |    lev0,lev1,lon0,lon1,lat)
!
! Save updated U and V:
        do i=lonbeg,lonend 
          do k=lev0,lev1-1
            un_upd(k,i,lat) = xx(1,k,i)
            vn_upd(k,i,lat) = xx(2,k,i)
!** No Wind **
!           un_upd(k,i,lat) = 0.
!           vn_upd(k,i,lat) = 0.
          enddo ! k=lev0,lev1-1
!
! Put spval in top nlevp1 level:
          un_upd(lev1,i,lat) = spval
          vn_upd(lev1,i,lat) = spval
!** No Wind **
!         un_upd(lev1,i,lat) = 0.
!         vn_upd(lev1,i,lat) = 0.
!
! Lower boundary is in ulbc,vlbc (fields.F):
          ulbc_nm(i,lat) = ulbc(i,lat)  ! LB for un_nm(itc) = un(itp)
          vlbc_nm(i,lat) = vlbc(i,lat)  ! LB for vn_nm(itc) = vn(itp)
          ulbc(i,lat) = unlbc(i)        ! LB for un(itc)
          vlbc(i,lat) = vnlbc(i)        ! LB for vn(itc)
!** No Wind **
!         ulbc_nm(i,lat) = 0.0         
!         vlbc_nm(i,lat) = 0.0 
!         ulbc(i,lat) =  0.0 
!         vlbc(i,lat) = 0.0 
        enddo ! i=lonbeg,lonend

!       call addfld('UN_SOLV',' ',' ',un_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('VN_SOLV',' ',' ',vn_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! End third latitude scan:
      enddo ! lat=lat0,lat1
!
! Set kut for wave filtering according to dlat (2.5 or 5.0):
!     call set_wave_filter(36,kut_5,nlat,kutt)
!
! Fourier smoothing of U and V:
!     call filter_uv(un_upd(lev0:lev1-1,:,:),lev0,lev1-1,lon0,lon1,
!    |  lat0,lat1,kutt)
!     call filter_uv(vn_upd(lev0:lev1-1,:,:),lev0,lev1-1,lon0,lon1,
!    |  lat0,lat1,kutt)
      call filter_uv(un_upd(lev0:lev1-1,:,:),lev0,lev1-1,lon0,lon1,
     |  lat0,lat1,kut)
      call filter_uv(vn_upd(lev0:lev1-1,:,:),lev0,lev1-1,lon0,lon1,
     |  lat0,lat1,kut)
!
! Fourth latitude scan:
      do lat=lat0,lat1

!       call addfld('UN_FILT',' ',' ',un_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('VN_FILT',' ',' ',vn_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Smooth updated un,vn:
        do i=lon0,lon1 
          do k=lev0,lev1-1
            unm_upd(k,i,lat) = dtsmooth_div2*(un_nm(k,i,lat)+
     |        un_upd(k,i,lat)) + dtsmooth*un(k,i,lat)
            vnm_upd(k,i,lat) = dtsmooth_div2*(vn_nm(k,i,lat)+
     |        vn_upd(k,i,lat)) + dtsmooth*vn(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1 

! Calculating Scale Height to convert winds from /s to m/s
        do i=lon0,lon1 
          do k=lev0+1,lev1-1
! hsc = scale height
           hsc(k,i) = gask*tni(k,i)/barm(k,i,lat)
           hsc(k,i) = hsc(k,i)/grav
          enddo ! k=lev0,lev1-1
           hsc(lev0,i) = gask*tni(lev0,i)/barm(lev0,i,lat)
           hsc(lev1,i) = gask*tni(lev1,i)/barm(lev1,i,lat)
           hsc(lev0,i) = hsc(lev0,i)/grav
           hsc(lev1,i) = hsc(lev1,i)/grav
        enddo ! i=lon0,lon1 
        do i=lon0,lon1 
          do k=lev0,lev1-1
           w_updm(k,i,lat) = w_upd(k,i,lat)*hsc(k,i)/100.
!         enddo ! k=lev0,lev1-1
!       enddo ! i=lon0,lon1 
!** Calculating Total drag term on the system tdrag * Utotal [m/sec2]**
!       do i=lon0,lon1
!         do k=lev0,lev1-1
            dragtermdcs(k,i,lat) = -(un_upd(k,i,lat)-
     |                              (2.0*cs(lat)*vrt2(k)))*tdrag(k)
            dragtermcs(k,i,lat) = -(un_upd(k,i,lat)-
     |                              (vrt2(k)*cs(lat)))*tdrag(k)
!
            dragtermdcs(k,i,lat) = dragtermdcs(k,i,lat)/100.
            dragtermcs(k,i,lat) = dragtermcs(k,i,lat)/100.
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!      do i = lonbeg,lonend
!        do k = lev0,lev1 - 1
!          write(6,*) ' un_upd( ',k,' ',i,' ) = ', un_upd(k,i,lat)
!          write(6,*) ' vn_upd( ',k,' ',i,' ) = ', vn_upd(k,i,lat)
!        enddo
!      enddo
!
!       call addfld('UN_NMOUT',' ',' ',unm_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('VN_NMOUT',' ',' ',vnm_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! End fourth latitude scan:
      enddo ! lat=lat0,lat1
!
#ifdef MPI
!
! Periodic points:
!     call mp_periodic_f3d(un_upd(:,lon0:lon1,lat0-1:lat1+1),
!    |  lev0,lev1,lon0,lon1,lat0-1,lat1+1)
!     call mp_periodic_f3d(vn_upd(:,lon0:lon1,lat0-1:lat1+1),
!    |  lev0,lev1,lon0,lon1,lat0-1,lat1+1)
#endif
      do lat=lat0,lat1
        call addfld('UN_FINAL',' ',' ',un_upd(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('VN_FINAL',' ',' ',vn_upd(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('WMSEC',' ',' ',w_updm(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('DRAG',' ',' ',dragtermcs(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('DRAGD',' ',' ',dragtermdcs(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      end subroutine duv
!-----------------------------------------------------------------------
      subroutine glp(tn,tn_upd,tn_nm,barm,z,zl,zp,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Horizontal pressure forcing for un,vn: 
!
      use params_module,only: nlonp4,nlat
      use cons_module,only: dz,dzgrav,freq_semidi,dt,grav,racs,re_inv
      use init_module,only: iter
!     use bndry_module,only: zb,zb2,zba,bnd,bnd2,bnda,ci
      use addfld_module,only: addfld
#ifdef MPI
      use mpi_module,only: mp_bndlons_f3d, mp_bndlats_f3d, 
     |  mp_periodic_f3d
#endif
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
!
! Inputs at full subdomain:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn,      ! neutral temperature (deg K)
     |  tn_upd,  ! updated neutral temperature (from dt.F)
     |  tn_nm,   ! neutral temperature from time n-1
     |  barm,    ! mean molecular weight
     |  z        ! geopotential
!
! Outputs at subdomain:
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(out) :: 
     |  zl,zp
!
! Local:
      integer :: k,i,nlevs,lat,lonbeg,lonend
      real,parameter :: wt=0.225
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2) ::
     |  ztmp
      real,dimension(lev0:lev1,lon0-2:lon1+2) ::
     |  tbar,barmi,dztbar
!
      nlevs = lev1-lev0+1
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = nlonp4-2
!
! Latitude scan:
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            tbar(k,i) = tn(k,i,lat)*(-2.)+tn_nm(k,i,lat)     ! s6
            tbar(k,i) = (tbar(k,i)+tn_upd(k,i,lat))*wt
            tbar(k,i) = tbar(k,i)+tn(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       if (lat <= lat1 .and. lat >= lat0) then
!         call addfld('TBAR',' ',' ',tbar(:,lon0:lon1),
!    |      'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       endif

        do i=lon0,lon1
          do k=lev0,lev1-1
            barmi(k,i) = .5*(barm(k,i,lat)+barm(k+1,i,lat))    ! s7
            dztbar(k,i) = (dz/dzgrav) * (tbar(k,i)/barmi(k,i)) ! s7
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       if (lat >= lat0 .and. lat <= lat1) then
!         call addfld('DZTBAR',' ',' ',dztbar(:,lon0:lon1),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       endif
!
! In earlier model versions, geopotential bottom boundary was recalculated
! here. In this code, it is just copied from z, where it is already available:
        do i=lon0,lon1
          ztmp(lev0,i,lat) = z(lev0,i,lat)
        enddo ! i=lon0,lon1
!
! This will differ from z(..[itp,itc]) because tn_upd is in dztbar.
        do i=lon0,lon1
          do k=lev0,lev1-1
            ztmp(k+1,i,lat) = dztbar(k,i)+ztmp(k,i,lat)
          enddo ! k=lev0,lev1-1

!         write(6,"('glp before bndlons: lat=',i2,' i=',i2,
!    |      ' lonbeg,end=',2i3,' ztmp(:,i,lat)=',/,(6e12.4))") 
!    |      lat,i,lonbeg,lonend,ztmp(:,i,lat)
        enddo ! i=lon0,lon1

!       if (lat >= lat0 .and. lat <= lat1) then
!         call addfld('ZTMP1',' ',' ',ztmp(:,lon0:lon1,lat),
!    |      'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       endif
      enddo ! lat=lat0,lat1
!
#ifdef MPI
!
! Exchange boundary longitudes in ztmp for 4th order differencing in dldp. 
! Not necessary if ztmp has been defined at full subdomain, including lon 
! halos, but as of 4/02, this was not working (s.a. comments above)
!
      call mp_bndlons_f3d(ztmp,nlevs,lon0,lon1,lat0-2,lat1+2,1)
      call mp_periodic_f3d(ztmp(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1)
#endif
!     do lat=lat0,lat1
!       call addfld('ZTMP2',' ',' ',ztmp(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)

!       do k=lev0,lev1
!         write(6,"('glp after bndlons: lat=',i2,' k=',i2,
!    |      ' ztmp(k,:,lat)=',/,(6e12.4))") 
!    |      lat,k,ztmp(k,:,lat)
!       enddo

!     enddo ! lat=lat0,lat1
!
! Calculate latitudinal and longitudinal derivatives in zl,zp:
! (dldp input ztmp is defined at full subdomain)
!
      call dldp(ztmp,z,zl,zp,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Complete zl,zp at subdomain:
!
      do lat=lat0,lat1
        do i=lonbeg,lonend
          do k=lev0,lev1-1
            zl(k,i,lat) = (zl(k,i,lat)+zl(k+1,i,lat))*
     |        (.5*grav*racs(lat))
            zp(k,i,lat) = (zp(k,i,lat)+zp(k+1,i,lat))*
     |        (.5*grav*re_inv)
          enddo ! k=lev0,lev1-1
        enddo ! i=lonbeg,lonend
!
! Set periodic points zero to avoid NaNS fpe in duv.
! (earlier versions actually set periodic points here)
        if (lon0==1) then
          zl(:,lon0:lon0+1,lat) = 0.
          zp(:,lon0:lon0+1,lat) = 0.
        endif
        if (lon1==nlonp4) then
          zl(:,lon1-1:lon1,lat) = 0.
          zp(:,lon1-1:lon1,lat) = 0.
        endif
      enddo ! lat=lat0,lat1

!     do lat=lat0,lat1
!       call addfld('ZL_GLP',' ',' ',zl(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('ZP_GLP',' ',' ',zp(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1

      end subroutine glp
!-----------------------------------------------------------------------
      subroutine dldp(zin,z,zl,zp,lev0,lev1,lon0,lon1,lat0,lat1)
      use params_module,only: nlonp4
      use cons_module,only: dlamda_1div12,dlamda_2div3,
     |  dphi_2div3,dphi_1div12
      use addfld_module,only: addfld
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
!
! Inputs at full subdomain:
! (Note new Z with new temperatures from glp is available only at the
!  current latitude, not adjacent latitudes. Therefore, the new Z 
!  from glp is used in longitudinal derivatives, whereas the current
!  Z (w/o updated temperature) is used in the latitudinal derivatives.
!  See ztmp in glp.)
!
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  zin,   ! new     Z (with updated temperature) from glp for zl
     |  z      ! current Z (with current temperature) from glp for zp
!
! Output at subdomain:
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(out) :: 
     |  zl,zp  ! output at subdomain
!
! Local:
      integer :: lonbeg,lonend,k,i,lat
!
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = nlonp4-2
!
! Longitudinal derivatives:
      do lat=lat0,lat1

!       call addfld('DLDP_ZIN',' ',' ',zin(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)

        do i=lonbeg,lonend
          do k=lev0,lev1
            zl(k,i,lat) = 
     |        dlamda_2div3 *(zin(k,i+1,lat)-zin(k,i-1,lat))-
     |        dlamda_1div12*(zin(k,i+2,lat)-zin(k,i-2,lat))
          enddo ! k=lev0,lev1
        enddo ! i=lonbeg,lonend
!
! Latitudinal derivatives:
        do i=lonbeg,lonend
          do k=lev0,lev1
            zp(k,i,lat) = 
     |        dphi_2div3 *(z(k,i,lat+1)-z(k,i,lat-1))-
     |        dphi_1div12*(z(k,i,lat+2)-z(k,i,lat-2))
          enddo ! k=lev0,lev1
        enddo ! i=lonbeg,lonend

!       call addfld('Z_JM2',' ',' ',z(:,lon0:lon1,lat-2),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('Z_JM1',' ',' ',z(:,lon0:lon1,lat-1),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('Z_JP1',' ',' ',z(:,lon0:lon1,lat+1),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('Z_JP2',' ',' ',z(:,lon0:lon1,lat+2),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)

!       call addfld('DLDP_ZL',' ',' ',zl(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('DLDP_ZP',' ',' ',zp(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)

      enddo ! lat=lat0,lat1
      end subroutine dldp
!-----------------------------------------------------------------------
      subroutine blktri(a,b,c,f,beta,gamma,y,x,lev0,lev1,lon0,lon1,lat)
!
! This procedure solves (I2-I1+1) tridiagonal block matrix
! systems in which all blocks are 2 x 2 matrices.
!
      use params_module,only: nlonp4 
      implicit none
!
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(2,2,lev0:lev1,lon0:lon1),intent(in) ::
     |  a,b,c
      real,dimension(2,lev0:lev1,lon0:lon1),intent(in) ::
     |  f
!
! Output args:
      real,dimension(2,2,lev0:lev1,lon0:lon1),intent(out) ::
     |  beta,gamma
      real,dimension(2,lev0:lev1,lon0:lon1),intent(out) ::
     |  x,y
!
! Local:
      integer :: k,i,lonbeg,lonend
!
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = nlonp4-2
!
! Lower boundary:
      do i=lonbeg,lonend
!
! Y(1,lev0,i) = determinant(B(lev0))
        y(1,lev0,i) = b(1,1,lev0,i)*b(2,2,lev0,i)-
     |                b(1,2,lev0,i)*b(2,1,lev0,i)
!
! BETA(K1) = B(K1)**(-1)
        beta(1,1,lev0,i) =  b(2,2,lev0,i)/y(1,lev0,i)
        beta(1,2,lev0,i) = -b(1,2,lev0,i)/y(1,lev0,i)
        beta(2,1,lev0,i) = -b(2,1,lev0,i)/y(1,lev0,i)
        beta(2,2,lev0,i) =  b(1,1,lev0,i)/y(1,lev0,i)
!
! Y(K1) = BETA(K1)*F(K1)
        y(1,lev0,i) = beta(1,1,lev0,i)*f(1,lev0,i)+
     |                beta(1,2,lev0,i)*f(2,lev0,i)
        y(2,lev0,i) = beta(2,1,lev0,i)*f(1,lev0,i)+
     |                beta(2,2,lev0,i)*f(2,lev0,i)
!
! Now deal with levels (K1+1),K2,1
        do k=lev0+1,lev1-1
!
! GAMMA(K-1) = BETA(K-1)*C(K-1)
          gamma(1,1,k-1,i) = beta(1,1,k-1,i)*c(1,1,k-1,i)+
     |                       beta(1,2,k-1,i)*c(2,1,k-1,i)          
          gamma(1,2,k-1,i) = beta(1,1,k-1,i)*c(1,2,k-1,i)+
     |                       beta(1,2,k-1,i)*c(2,2,k-1,i)          
          gamma(2,1,k-1,i) = beta(2,1,k-1,i)*c(1,1,k-1,i)+
     |                       beta(2,2,k-1,i)*c(2,1,k-1,i)          
          gamma(2,2,k-1,i) = beta(2,1,k-1,i)*c(1,2,k-1,i)+
     |                       beta(2,2,k-1,i)*c(2,2,k-1,i)          
!
! GAMMA(K) = B(K) - A(K)*GAMMA(K-1)
          gamma(1,1,k,i) = b(1,1,k,i)-a(1,1,k,i)*gamma(1,1,k-1,i)-
     |                                a(1,2,k,i)*gamma(2,1,k-1,i)
          gamma(1,2,k,i) = b(1,2,k,i)-a(1,1,k,i)*gamma(1,2,k-1,i)-
     |                                a(1,2,k,i)*gamma(2,2,k-1,i)
          gamma(2,1,k,i) = b(2,1,k,i)-a(2,1,k,i)*gamma(1,1,k-1,i)-
     |                                a(2,2,k,i)*gamma(2,1,k-1,i)
          gamma(2,2,k,i) = b(2,2,k,i)-a(2,1,k,i)*gamma(1,2,k-1,i)-
     |                                a(2,2,k,i)*gamma(2,2,k-1,i)
!
! Y(1,k,i) = determinant(GAMMA(K))
          y(1,k,i) = gamma(1,1,k,i)*gamma(2,2,k,i)-
     |               gamma(1,2,k,i)*gamma(2,1,k,i)
!
! BETA(K) = GAMMA(K)**(-1)
          beta(1,1,k,i) =  gamma(2,2,k,i)/y(1,k,i)
          beta(1,2,k,i) = -gamma(1,2,k,i)/y(1,k,i)
          beta(2,1,k,i) = -gamma(2,1,k,i)/y(1,k,i)
          beta(2,2,k,i) =  gamma(1,1,k,i)/y(1,k,i)
!
!  X(K) = F(K) - A(K)*Y(K-1)
          x(1,k,i) = f(1,k,i)-a(1,1,k,i)*y(1,k-1,i)-
     |                        a(1,2,k,i)*y(2,k-1,i)
          x(2,k,i) = f(2,k,i)-a(2,1,k,i)*y(1,k-1,i)-
     |                        a(2,2,k,i)*y(2,k-1,i)
!
! Y(K) = BETA(K)*X(K)
          y(1,k,i) = beta(1,1,k,i)*x(1,k,i)+beta(1,2,k,i)*x(2,k,i)
          y(2,k,i) = beta(2,1,k,i)*x(1,k,i)+beta(2,2,k,i)*x(2,k,i)
        enddo ! k=lev0+1,lev1-1
      enddo ! i=lonbeg,lonend
!
      do i=lonbeg,lonend
!
!  X(K2) = Y(K2)
        x(1,lev1-1,i) = y(1,lev1-1,i)
        x(2,lev1-1,i) = y(2,lev1-1,i)
!
! Backward sweep to determine final solution, X(K) for k=k2,k1,-1
        do k=lev1-2,lev0,-1
          x(1,k,i) = y(1,k,i)-gamma(1,1,k,i)*x(1,k+1,i)-
     |                        gamma(1,2,k,i)*x(2,k+1,i)
          x(2,k,i) = y(2,k,i)-gamma(2,1,k,i)*x(1,k+1,i)-
     |                        gamma(2,2,k,i)*x(2,k+1,i)
        enddo ! k=lev1-2,lev1-1,-1
      enddo ! i=lon0,lon1
      end subroutine blktri
!-----------------------------------------------------------------------
      subroutine filter_uv(fout,lev0,lev1,lon0,lon1,lat0,lat1,kut)
!
! Filter updated W omega:
!
      use params_module,only: nlat,nlonp4,nlon
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
!     code = 131 ; state = 'filter_uv' ; activity='Filtering'
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
!     code = 131 ; state = 'filter_uv' ; activity='Filtering'
      call vtend(131,ier)
#endif
      end subroutine filter_uv
!-----------------------------------------------------------------------
      subroutine calc_uvzm(un,uzm,vn,vzm,lev0,lev1,lon0,lon1,lat0,lat1)
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
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(in) :: un,vn
      real,dimension(lev0:lev1,lat0:lat1),intent(out) :: uzm,vzm
!
! Local:
      integer :: k,i,j
      real :: f_kij(lev0:lev1,nlonp4,lat0:lat1,2)
      real :: rlon

      f_kij = 0.
      do j=lat0,lat1
        do i=lon0,lon1
          f_kij(:,i,j,1) = un(:,i,j)
          f_kij(:,i,j,2) = vn(:,i,j)
        enddo
      enddo

#ifdef MPI
!
! Gather global longitudes to western-most task in lat row:
      call mp_gatherlons_f3d(f_kij,lev0,lev1,lon0,lon1,lat0,lat1,2)
#endif
      uzm = 0. ; vzm = 0.
      rlon = 1./real(nlonp4-4)
!
! Calculate zonal mean:
      if (mytidi==0) then
        do j=lat0,lat1
          do i=3,nlonp4-2
            uzm(:,j) = uzm(:,j)+f_kij(:,i,j,1)
            vzm(:,j) = vzm(:,j)+f_kij(:,i,j,2)
          enddo
          uzm(:,j) = uzm(:,j)*rlon
          vzm(:,j) = vzm(:,j)*rlon
        enddo
!
! Copy zm redundantly to i-index of f_kij:
        do j=lat0,lat1
          do i=1,nlonp4
            f_kij(:,i,j,1) = uzm(:,j)
            f_kij(:,i,j,2) = vzm(:,j)
          enddo
        enddo
      endif ! mytidi==0
#ifdef MPI
!
! Scatter 3d w/ zonal means to other tasks in my lat row:
      call mp_scatterlons_f3d(f_kij,lev0,lev1,lon0,lon1,lat0,lat1,2)
#endif
!
! Finally copy zm to my subdomain tzm:
      do j=lat0,lat1
        uzm(:,j) = f_kij(:,lon0,j,1) ! all i's in f_kij are zm, so just use lon0
        vzm(:,j) = f_kij(:,lon0,j,2)
      enddo
      end subroutine calc_uvzm
!-----------------------------------------------------------------------
