!
      subroutine comp(tn,difk,o1,o1_nm,co,co_nm,n2,un,vn,w,hdo1,hdco,
     |                o1_upd,o1nm_upd,co_upd,conm_upd,
     |                lev0,lev1,lon0,lon1,lat0,lat1)
!
! Advance major species O and CO.
! 10/06 btf: Adapted for vtgcm, from tiegcm1.8 
!            (replace o2 with co, and o1 with o2)
! 11/17/06 swb: corrections for vtgcm physics (n2 terms added)
! 12/04/06 swb: added difk to arg list (from edddyflds)
! 01/30/08 swb: recheck ps0(:,1 and 2), embar0(:) and wks4(:)
! 05/14/08 swb: small = 1.0E-10 to enable PSI01 to reach specified lbc.
! 05/15/08 swb: small = 1.0E-11 to enable PSI01 to reach specified lbc.
! 05/15/08 swb: small = 1.0E-13 to enable PSI01 to reach specified lbc.
!
      use params_module,only: dz,nlonp4,nlat,dlat,nlevp1
      use input_module,only: difhor
      use init_module,only: glat,istep
      use cons_module,only: pi,hor,dtr,rmassinv_o1,rmassinv_co,
     |  rmassinv_n2,rmass_o1,rmass_co,expz,expzmid,expzmid_inv,
     |  dtx2inv,kut_5,kut,dtsmooth,dtsmooth_div2,
     |  rmassinv_co2,rmass_n2,rmass_co2
      use chemrates_module,only: fs
      use lbc,only: b,fb
      use addfld_module,only: addfld
      use chapman_module,only: idn_lonlat
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn,   ! neutral temperature
     |  difk, ! eddy diffusion coefficient (from eddyflds)
     |  o1,   ! O  (mmr) at current timestep
     |  co,   ! CO  (mmr) at current timestep
     |  n2,   ! N2  (mmr) at current timestep
     |  o1_nm,! O  (mmr) at time n-1
     |  co_nm,! CO  (mmr) at time n-1
     |  un,   ! zonal wind velocity at current timestep
     |  vn,   ! meridional wind velocity at current timestep
     |  w,    ! vertical velocity at current timestep
     |  hdo1, ! O  horizontal diffusion (hdif3)
     |  hdco  ! CO  horizontal diffusion (hdif3)
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out)::
     |  o1_upd   ,co_upd,  ! output: O,CO updated for next timestep
     |  o1nm_upd ,conm_upd ! output: O,CO updated for previous timestep
!
! Local:
      integer :: k,kk,i,lat,isp,km,kp,ktmp,m,lonbeg,lonend,kutt(nlat)
      integer :: nk,nkm1
      integer,parameter :: io1=1,ico=2 ! indices to O, CO, respectively
      real,dimension(lon0:lon1,lev0:lev1,2,2) :: gama
      real,dimension(lon0:lon1,lev0:lev1,2)   :: zz
      real,dimension(lon0:lon1,lev0:lev1)     :: embar
      real,dimension(lon0:lon1,2,2,2)         :: ak
      real,dimension(lon0:lon1,2,2)           :: ep,pk,qk,rk,wkm1,wkm2
      real,dimension(lon0:lon1,2)             :: fk,wkv1,wkv2,ps0
      real,dimension(lon0:lon1)               :: wks1,wks2,wks3,wks4,
     |                                           embar0,dfactor
      real :: rlat
      real :: ak0(2,2),phi(2,3),delta(2,2),tau,t00,small
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2)::
     |  oref,co2ref,coref,n2ref                                ! Ref. mmr for testing
      real,dimension(nlevp1)  :: orefd,co2refd,corefd,n2refd,  ! Ref. vmr for testing
     |  orefn,co2refn,corefn,n2refn,mbarrefd,mbarrefn                         
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  o1nm_smooth, conm_smooth,  ! smoothed at time n-1
     |  o1_advec   , co_advec,     ! horizontal advection
     |  psco2                      ! YCO2 (local)
      real,dimension(lev0:lev1,lon0:lon1,2) :: upd
      logical,parameter :: debug=.false.
!
! For diagnostics:
      real,dimension(lev0:lev1,lon0:lon1,2) :: zz_ki
      real,dimension(lev0:lev1,lon0:lon1,2,2) :: gama_ki

      if (debug) write(6,"('Enter comp: istep=',i3)") istep
!
! Earth:
!     phi(:,1)=(/0.  ,0.673/)
!     phi(:,2)=(/1.35,0.   /)
!     phi(:,3)=(/1.11,0.769/)
!     tau=1.86e+3
!     delta(:,1)=(/1.,0./)
!     delta(:,2)=(/0.,1./)
!     t00=273.
!
! Venus:
      phi(:,1) = (/0.    ,2.296/)
      phi(:,2) = (/1.31  ,0.   /)
      phi(:,3) = (/1.026 ,1.600/)
      tau=9.25E+3
      delta(:,1)=(/1., 0./)
      delta(:,2)=(/0., 1./)
      t00=273.
!
!     small = 1.e-6
!     small = 1.e-11
      small = 1.e-13
      nk = lev1-lev0+1
      nkm1 = nk-1
!
! Calculate and save horizontal advection in o1_advec, co_advec:
!
      do lat=lat0,lat1
        call advecl(o1,co,un,vn,o1_advec,co_advec,
     |    lev0,lev1,lon0,lon1,lat0,lat1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('comp after advecl: istep=',i3)") istep
!
! Save smoothed o1,o at time n-1:
!
      call smooth(o1_nm,o1nm_smooth,lev0,lev1,lon0,lon1,lat0,lat1,0)
      if (debug) write(6,"('comp after smooth O: istep=',i3)") istep
      call smooth(co_nm,conm_smooth,lev0,lev1,lon0,lon1,lat0,lat1,0)
      if (debug) write(6,"('comp after smooth CO: istep=',i3)") istep

      do lat=lat0,lat1
        call addfld('O1SMOOTH',' ',' ',o1nm_smooth(:,:,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('COSMOOTH',' ',' ',conm_smooth(:,:,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
! Begin latitude scan:
      do lat=lat0,lat1
!
!****Miguel Reference Atmosphere Data Blocks****
! Composition is VMR and then converted to MMR

! ** DAY **

!     orefd(:) = (/
!    | 2.7374e-09,5.4469e-09,7.5505e-09,1.0504e-08,1.3992e-08,
!    | 1.1221e-08,7.5341e-09,1.4539e-08,5.2946e-08,1.9000e-07,
!    | 6.6008e-07,2.1737e-06,7.0352e-06,1.8487e-05,3.5288e-05,
!    | 6.4598e-05,1.0666e-04,1.5665e-04,2.2818e-04,3.3720e-04,
!    | 5.0233e-04,7.5314e-04,1.1479e-03,1.7397e-03,2.6943e-03,
!    | 4.1820e-03,6.5031e-03,1.0131e-02,1.5670e-02,2.4151e-02,
!    | 3.7113e-02,5.6169e-02,8.3620e-02,1.2150e-01,1.7075e-01,
!    | 2.3185e-01,3.0356e-01,3.8313e-01,4.6700e-01,5.5188e-01,
!    | 6.3405e-01,7.0965e-01,7.7452e-01,8.2799e-01,
!    | 8.2799e-01,8.2799e-01,8.2799e-01,8.2799e-01,8.2799e-01,
!    | 8.2799e-01,8.2799e-01,8.2799e-01,8.2799e-01,8.2799e-01,
!    | 8.2799e-01,8.2799e-01,8.2799e-01,8.2799e-01,8.2799e-01,
!    | 8.2799e-01,8.2799e-01,8.2799e-01,8.2799e-01,8.2799e-01,
!    | 8.2799e-01,8.2799e-01,8.2799e-01,8.2799e-01,8.2799e-01/)
!
!     co2refd(:) = (/
!    | 0.9650,0.9650,0.9650,0.9650,0.9650,0.9650,0.9650,0.9650,
!    | 0.9650,0.9650,0.9632,0.9624,0.9615,0.9607,0.9598,0.9590,
!    | 0.9580,0.9571,0.9570,0.9543,0.9513,0.9500,0.9491,0.9490,
!    | 0.9416,0.9341,0.9236,0.9113,0.8952,0.8728,0.8460,0.8052,
!    | 0.7561,0.6931,0.6144,0.5299,0.4407,0.3496,0.2654,0.1885,
!    | 0.1269,0.0784,0.0449,0.0245,
!    | 0.0245,0.0245,0.0254,0.0254,0.0254,0.0254,
!    | 0.0245,0.0245,0.0254,0.0254,0.0254,0.0254,
!    | 0.0245,0.0245,0.0254,0.0254,0.0254,0.0254,
!    | 0.0245,0.0245,0.0254,0.0254,0.0254,0.0254,
!    | 0.0245/)
!
!     corefd(:) = (/
!    | 2.0143e-06,2.0000e-06,2.0000e-06,2.0000e-06,2.0000e-06,
!    | 2.0000e-06,2.0000e-06,2.0000e-06,2.1145e-06,3.0311e-06,
!    | 5.4432e-06,1.0167e-05,1.8829e-05,3.4894e-05,6.4560e-05,
!    | 1.1368e-04,1.8300e-04,2.7225e-04,3.9878e-04,5.8882e-04,
!    | 8.7047e-04,1.2928e-03,1.9209e-03,2.8596e-03,4.2632e-03,
!    | 6.3543e-03,9.4045e-03,1.3911e-02,2.0221e-02,2.9202e-02,
!    | 4.1307e-02,5.7700e-02,7.7858e-02,1.0184e-01,1.2839e-01,
!    | 1.5407e-01,1.7514e-01,1.8971e-01,1.9500e-01,1.9062e-01,
!    | 1.7784e-01,1.5674e-01,1.3156e-01,1.0563e-01,
!    | 1.0563e-01,1.0563e-01,1.0563e-01,1.0563e-01,1.0563e-01,
!    | 1.0563e-01,1.0563e-01,1.0563e-01,1.0563e-01,1.0563e-01,
!    | 1.0563e-01,1.0563e-01,1.0563e-01,1.0563e-01,1.0563e-01,
!    | 1.0563e-01,1.0563e-01,1.0563e-01,1.0563e-01,1.0563e-01,
!    | 1.0563e-01,1.0563e-01,1.0563e-01,1.0563e-01,1.0563e-01/)
!
!      mbarrefd(:) = (/
!    | 44.,44.,44.,44.,44.,44.,44.,44.,44.,44.,
!    | 43.9,43.9,43.9,43.9,43.8,43.5,43.3,43.3,
!    | 43.3,43.3,43.2,43.2,43.2,43.1,43.0,42.9,
!    | 42.7,42.5,42.2,41.8,41.2,40.3,39.2,37.7,
!    | 35.9,33.7,31.4,29.0,26.6,24.4,22.4,20.6,
!    | 20.6,20.6,20.6,20.6,20.6,20.6,20.6,20.6,
!    | 20.6,20.6,20.6,20.6,20.6,20.6,20.6,20.6,
!    | 20.6,20.6,20.6,20.6,20.6,20.6,20.6,20.6,
!    | 20.6,20.6,20.6/)
!
!
! ** NIGHT **
!
!      orefn(:) = (/
!    | 2.7374e-09,5.4469e-09,7.5505e-09,1.0504e-08,1.3943e-08,
!    | 1.0368e-08,6.6430e-09,1.6455e-08,5.7191e-08,1.9586e-07,
!    | 6.5912e-07,2.1959e-06,7.3556e-06,1.7827e-05,3.3027e-05,
!    | 5.9893e-05,9.6665e-05,1.3532e-04,1.8710e-04,2.6334e-04,
!    | 3.7051e-04,5.2018e-04,7.4533e-04,1.0650e-03,1.5145e-03,
!    | 2.2032e-03,3.1767e-03,4.6325e-03,6.7085e-03,9.7297e-03,
!    | 1.4145e-02,2.0446e-02,2.9672e-02,4.2795e-02,6.1617e-02,
!    | 8.7124e-02,1.2390e-01,1.7312e-01,2.3426e-01,3.1736e-01,
!    | 4.1816e-01,5.2107e-01,6.4465e-01,7.3253e-01,8.1405e-01,
!    | 8.6450e-01,8.8982e-01,8.9957e-01,8.9290e-01,8.6728e-01,
!    | 8.2019e-01,7.5296e-01,6.6145e-01,5.3957e-01,3.9177e-01,
!    | 3.9177e-01,3.9177e-01,3.9177e-01,3.9177e-01,3.9177e-01,
!    | 3.9177e-01,3.9177e-01,3.9177e-01,3.9177e-01,3.9177e-01,
!    | 3.9177e-01,3.9177e-01,3.9177e-01,3.9177e-01/)
!
!      co2refn(:) = (/
!    | 0.9650,0.9650,0.9650,0.9650,0.9650,0.9650,0.9650,0.9650,
!    | 0.9650,0.9650,0.9650,0.9650,0.9650,0.9650,0.9650,0.9650,
!    | 0.9635,0.9618,0.9600,0.9583,0.9567,0.9550,0.9527,0.9504,
!    | 0.9481,0.9450,0.9418,0.9387,0.9296,0.9188,0.9092,0.8921,
!    | 0.8706,0.8492,0.8206,0.7694,0.7216,0.6730,0.5705,0.4828,
!    | 0.3894,0.2717,0.1873,0.1051,0.0565,2.7268e-02,
!    | 1.2307e-02,5.2888e-03,2.1558e-03,8.4015e-04,
!    | 3.0147e-04,1.0085e-04,2.9519e-05,7.0140e-06,1.2316e-06,
!    | 1.2316e-06,1.2316e-06,1.2316e-06,1.2316e-06,1.2316e-06,
!    | 1.2316e-06,1.2316e-06,1.2316e-06,1.2316e-06,1.2316e-06,
!    | 1.2316e-06,1.2316e-06,1.2316e-06,1.2316e-06/)
!
!      corefn(:) = (/
!    | 2.0143e-06,2.0000e-06,2.0000e-06,2.0000e-06,2.0000e-06,
!    | 2.0000e-06,2.0000e-06,2.0000e-06,2.0954e-06,2.8335e-06,
!    | 4.7427e-06,8.2258e-06,1.4380e-05,2.4843e-05,4.3164e-05,
!    | 7.0236e-05,1.0077e-04,1.3832e-04,1.8912e-04,2.6072e-04,
!    | 3.5913e-04,4.9483e-04,6.8273e-04,9.4088e-04,1.2996e-03,
!    | 1.7887e-03,2.4634e-03,3.3996e-03,4.6208e-03,6.2741e-03,
!    | 8.4565e-03,1.1287e-02,1.5017e-02,1.9777e-02,2.5673e-02,
!    | 3.2886e-02,4.1384e-02,5.0643e-02,6.0534e-02,7.0026e-02,
!    | 7.7357e-02,8.0528e-02,7.9128e-02,7.0728e-02,5.9386e-02,
!    | 4.6159e-02,3.4054e-02,2.4220e-02,1.6519e-02,1.0921e-02,
!    | 6.8505e-03,4.0950e-03,2.2527e-03,1.0844e-03,4.3112e-04,
!    | 4.3112e-04,4.3112e-04,4.3112e-04,4.3112e-04,4.3112e-04,
!    | 4.3112e-04,4.3112e-04,4.3112e-04,4.3112e-04,4.3112e-04,
!    | 4.3112e-04,4.3112e-04,4.3112e-04,4.3112e-04/)
!
!      mbarrefn(:) = (/
!    | 44.,44.,44.,44.,44.,44.,44.,44.,44.,44.,
!    | 43.9,43.9,43.8,43.7,43.6,43.3,43.1,43.2,
!    | 43.2,43.2,43.1,43.1,43.2,43.2,43.0,42.7,
!    | 42.6,42.4,42.2,41.8,41.4,40.8,40.0,38.5,
!    | 37.1,34.5,32.0,29.7,27.7,26.2,25.7,26.6,
!    | 26.6,26.6,26.6,26.6,26.6,26.6,26.6,26.6,
!    | 26.6,26.6,26.6,26.6,26.6,26.6,26.6,26.6,
!    | 26.6,26.6,26.6,26.6,26.6,26.6,26.6,26.6,
!    | 26.6,26.6,26.6/)
!
!** convert vmr to mmr

!       do i=lon0,lon1   !1,nlonp4
!         do k=lev0,lev1   !1,nlevp1
!** Day-Night Variations Specified
!           if (idn_lonlat(i,lat)==1) then
!             oref(k,i,lat) = orefd(k)*rmass_o1/mbarrefd(k)
!             coref(k,i,lat) = corefd(k)*rmass_co/mbarrefd(k)
!             co2ref(k,i,lat) = co2refd(k)*rmass_co2/mbarrefd(k)
!           else
!             oref(k,i,lat) = orefn(k)*rmass_o1/mbarrefn(k)
!             coref(k,i,lat) = corefn(k)*rmass_co/mbarrefn(k)
!             co2ref(k,i,lat) = co2refn(k)*rmass_co2/mbarrefn(k)
!           endif
!         enddo
!       enddo


      do i=lon0,lon1
        dfactor(i) = 1.
      enddo
!
! If difhor flag was by user (see input_module):
! dfactor = 
!   .5*(1.+SIN(PI*(ABS(RLATM)-PI/6.)/(PI/3.))) FOR ABS(RLATM).LT.PI/3.
! dfactor = 1. FOR  ABS(RLATM).GE.PI/3
! (dfactor was in sub dfact in earlier versions)
!
      if (difhor > 0) then
        rlat = (glat(1)+(lat-1)*dlat)*dtr
        if (abs(rlat)-pi/4.5 >= 0.) then
          dfactor(:) = hor(lat)+1.
        else
          dfactor(:) = hor(lat)+.5*(1.+sin(pi*(abs(rlat)-pi/9.)/
     |      (pi/4.5)))
        endif
      else
        dfactor(:) = 1.
      endif
!     write(6,"('comp: lat=',i2,' dfactor=',/,(6e12.4))")
!    |  lat,dfactor
!
! Embar: Corrected for Venus: swb  (11/17/06)
      do i=lon0,lon1
        do k=lev0,lev1
	  psco2(k,i,lat)=  (1.-o1(k,i,lat)-co(k,i,lat)-n2(k,i,lat))
	  if (psco2(k,i,lat) < 0.) psco2(k,i,lat) = 0.
          embar(i,k) = 1./(o1(k,i,lat)*rmassinv_o1 + 
     |                     co(k,i,lat)*rmassinv_co +
     |                     n2(k,i,lat)*rmassinv_n2 +
     |                     psco2(k,i,lat)*rmassinv_co2)
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1
!
! Note embar is (i,k), not (k,i):
      call addfld('EMBAR',' ',' ',embar,
     |    'lon',lon0,lon1,'lev',lev0,lev1,lat)
!
! ps0 and embar0:
! Correct to include N2: swb (11/17/06)
! Recheck for accuracy (prints) to confirm proper calculation (1/30/08)
      do i=lon0,lon1
        ps0(i,io1) = 
     |    b(i,1,1)*o1(lev0,i,lat)+b(i,1,2)*co(lev0,i,lat)+fb(i,1) 
        ps0(i,ico) = 
     |    b(i,2,1)*o1(lev0,i,lat)+b(i,2,2)*co(lev0,i,lat)+fb(i,2) 
!       embar0(i) = 1./(ps0(i,io1)*rmassinv_o1+ps0(i,ico)*rmassinv_co+
!    |              (1.-ps0(i,io1)-ps0(i,ico))*rmassinv_n2)
        embar0(i) = 1./(ps0(i,io1)*rmassinv_o1+ps0(i,ico)*rmassinv_co+
     |        n2(lev0,i,lat)*rmassinv_n2+
     |        (1.-ps0(i,io1)-ps0(i,ico)-n2(lev0,i,lat))*rmassinv_co2)
! WKS4 = .5*(DMBAR/DZ)/MBAR
        wks4(i) = (embar(i,lev0)-embar0(i))/
     |            (dz*(embar0(i)+embar(i,lev0)))
      enddo ! i=lon0,lon1

! Recheck for accuracy (prints) to confirm proper calculation (1/30/08)
      if (debug) then
        write(6,"(/,'comp: lat=',i2,' ps0(:,io1)=',/,(6e12.4))") 
     |    lat,ps0(:,io1)
        write(6,"('ps0(:,ico)=',/,(6e12.4))") ps0(:,ico)
        write(6,"('embar0(:)=',/,(6e12.4))") embar0(:)
        write(6,"('wks4(:)=',/,(6e12.4))") wks4(:)
      endif
!
! ep, ak at level 1/2:
      km = 1 
      kp = 2
      do i=lon0,lon1
        ep(i,io1,kp) = 1.-(2./(embar0(i)+embar(i,lev0)))*
     |               (rmass_o1+(embar(i,lev0)-embar0(i))/dz)
        ep(i,ico,kp) = 1.-(2./(embar0(i)+embar(i,lev0)))*
     |               (rmass_co+(embar(i,lev0)-embar0(i))/dz)
        zz(i,1,:) = 0.
      enddo ! i=lon0,lon1

      if (debug) then
        write(6,"('comp: lat=',i3,' ep(:,1,kp)=',/,(6e12.4))") 
     |    lat,ep(:,1,kp)
        write(6,"('comp: lat=',i3,' ep(:,2,kp)=',/,(6e12.4))") 
     |    lat,ep(:,2,kp)
      endif

      do m=1,2
        do i=lon0,lon1
          ak(i,io1,m,kp) = 
     |      -delta(io1,m)*(phi(ico,3)+(phi(ico,io1)-phi(ico,3))*
     |      .5*(ps0(i,io1)+o1(lev0,i,lat)))-(1.-delta(io1,m))*
     |      (phi(io1,m)-phi(io1,3))*.5*(ps0(i,io1)+o1(lev0,i,lat))
          ak(i,ico,m,kp) = 
     |      -delta(ico,m)*(phi(io1,3)+(phi(io1,ico)-phi(io1,3))*
     |      .5*(ps0(i,ico)+co(lev0,i,lat)))-(1.-delta(ico,m))*
     |      (phi(ico,m)-phi(ico,3))*.5*(ps0(i,ico)+co(lev0,i,lat))
        enddo ! i=lon0,lon1

!       write(6,"('lat=',i2,' m=',i2,' ak(io1)=',/,(6e12.4))")
!    |    lat,m,ak(:,io1,m,kp)
!       write(6,"('lat=',i2,' m=',i2,' ak(ico)=',/,(6e12.4))")
!    |    lat,m,ak(:,ico,m,kp)

      enddo ! m=1,2
!
! WKS1=MBAR/M3*(T00/(T0+T))*0.25/(TAU*DET) AT LEVEL 1/2
! tn lower boundary is stored in top slot tn(lev1..). 
! Replace rmassinv_n2 with rmassinv_co2
      do i=lon0,lon1
        wks1(i) = 0.5*(embar0(i)+embar(i,lev0))*rmassinv_co2*
     |    (t00/tn(lev1,i,lat))**0.25/(tau*(ak(i,1,1,kp)*ak(i,2,2,kp)-
     |    ak(i,1,2,kp)*ak(i,2,1,kp)))
      enddo ! i=lon0,lon1
!     write(6,"('comp: lat=',i3,' wks1=',/,(6e12.4))") lat,wks1
!
! Complete claculation of ak(1/2)
      do m=1,2
        do i=lon0,lon1
          ak(i,io1,m,kp) = ak(i,io1,m,kp)*wks1(i)
          ak(i,ico,m,kp) = ak(i,ico,m,kp)*wks1(i)
          gama(i,lev0,:,m) = 0.
        enddo ! i=lon0,lon1

!       write(6,"('lat=',i2,' m=',i2,' ak(io1)=',/,(6e12.4))")
!    |    lat,m,ak(:,io1,m,kp)
!       write(6,"('lat=',i2,' m=',i2,' ak(ico)=',/,(6e12.4))")
!    |    lat,m,ak(:,ico,m,kp)

      enddo ! m=1,2
!
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = nlonp4-2
!
! Height (pressure) loop:
! For now (4/02), put k-loop on outside even tho embar and input
!   fields are (k,i), for convenience in verification with tgcm15.
!
      km = 1  ! alternates 2,1,2,1,... during k-loop
      kp = 2  ! alternates 1,2,1,2,... during k-loop
      levloop: do k=lev0,lev1-1                               ! DO 6
        ktmp = km
        km   = kp
        kp   = ktmp
        do i=lon0,lon1
          ep(i,io1,kp) = 1.-(2./(embar(i,k)+embar(i,k+1)))*(rmass_o1+
     |      (embar(i,k+1)-embar(i,k))/dz)
          ep(i,ico,kp) = 1.-(2./(embar(i,k)+embar(i,k+1)))*(rmass_co+
     |      (embar(i,k+1)-embar(i,k))/dz)
        enddo ! i=lon0,lon1

!       write(6,"('lat=',i3,' k=',i3,' ep(io1)=',/,(6e12.4))") lat,k,
!    |    ep(:,io1,kp)
!       write(6,"('ep(ico)=',/,(6e12.4))") ep(:,ico,kp)

        do m=1,2
          do i=lon0,lon1
!
! AK(K+1/2)
            ak(i,io1,m,kp) = 
     |        -delta(io1,m)*(phi(ico,3)+(phi(ico,io1)-phi(ico,3))*
     |        .5*(o1(k,i,lat)+o1(k+1,i,lat)))-
     |        (1.-delta(io1,m))*(phi(io1,m)-phi(io1,3))*
     |        .5*(o1(k,i,lat)+o1(k+1,i,lat))

            ak(i,ico,m,kp) = 
     |        -delta(ico,m)*(phi(io1,3)+(phi(io1,ico)-phi(io1,3))*
     |        .5*(co(k,i,lat)+co(k+1,i,lat)))-
     |        (1.-delta(ico,m))*(phi(ico,m)-phi(ico,3))*
     |        .5*(co(k,i,lat)+co(k+1,i,lat))

          enddo ! i=lon0,lon1
        enddo ! m=1,2
!
! WKS1=MBAR/M3*(T00/(T0+T))**0.25/(TAU*DET(ALFA))
! Replace rmassinv_n2 with rmassinv_co2
        do i=lon0,lon1
          wks1(i) = 0.5*(embar(i,k)+embar(i,k+1))*rmassinv_co2*
     |      (t00/(.5*(tn(k,i,lat)+tn(k+1,i,lat))))**0.25/
     |      (tau*(ak(i,1,1,kp)*ak(i,2,2,kp)-ak(i,1,2,kp)*
     |      ak(i,2,1,kp)))
!
! EDDY DIFFUSION TERMS IN WKS3 AND WKS4 
          wks3(i) = wks4(i)
          wks4(i) = (embar(i,k+1)-embar(i,k))/
     |      (dz*(embar(i,k)+embar(i,k+1)))
        enddo ! i=lon0,lon1
!
! FINISH CALCULATING AK(K+1/2) AND GENERATE PK, QK, RK
! DIFK (1/sec) used for eddy diffusion (see eddyflds.F)
        do m=1,2
          do isp=io1,ico
            do i=lon0,lon1
              ak(i,isp,m,kp) = ak(i,isp,m,kp)*wks1(i)

              pk(i,isp,m) = (ak(i,isp,m,km)*(1./dz+ep(i,m,km)/2.)-
     |          expz(k)*(expzmid_inv*difk(k,i,lat)*dfactor(i)*(1./dz-
     |          wks3(i))+0.25*(w(k,i,lat)+w(k+1,i,lat)))*
     |          delta(isp,m))/dz

              rk(i,isp,m) = (ak(i,isp,m,kp)*(1./dz-ep(i,m,kp)/2.)-
     |          expz(k)*(expzmid*difk(k+1,i,lat)*dfactor(i)*(1./dz+
     |          wks4(i))-0.25*(w(k,i,lat)+w(k+1,i,lat)))*
     |          delta(isp,m))/dz
   
              qk(i,isp,m) = -(ak(i,isp,m,km)*(1./dz-ep(i,m,km)/2.)+
     |          ak(i,isp,m,kp)*(1./dz+ep(i,m,kp)/2.))/dz+expz(k)*
     |          (((expzmid*difk(k+1,i,lat)*(1./dz-wks4(i))+expzmid_inv*
     |          difk(k,i,lat)*(1./dz+wks3(i)))*dfactor(i)/dz+dtx2inv)*
     |          delta(isp,m)-fs(i,k,isp,m,lat))

            enddo ! i=lon0,lon1

!           write(6,"(/,'comp: lat=',i3,' m=',i2,' isp=',i2,' k=',i3)") 
!    |        lat,m,isp,k
!           write(6,"('ak=',/,(6e12.4))") ak(:,isp,m,kp)
!           write(6,"('pk=',/,(6e12.4))") pk(:,isp,m)
!           write(6,"('rk=',/,(6e12.4))") rk(:,isp,m)
!           write(6,"('qk=',/,(6e12.4))") qk(:,isp,m)

          enddo ! isp=io1,ico
        enddo ! m=1,2
!
! Use advection saved from advecl calls at beginning of routine:
!
        do i=lon0,lon1
          fk(i,io1) = o1_advec(k,i,lat)
          fk(i,ico) = co_advec(k,i,lat)
        enddo ! i=lonbeg,lonend

!       write(6,"(/,'comp: lat=',i3,' k=',i2)") lat,k
!       write(6,"('fk(:,io1)=',/,(6e12.4))") fk(lonbeg:lonend,io1)
!       write(6,"('fk(:,ico)=',/,(6e12.4))") fk(lonbeg:lonend,ico)
!
! Add explicit source terms to fk: from comp_oco.F
! Turn on or off for testing:  November 2006
        do i=lon0,lon1
          fk(i,io1) = fk(i,io1)-fs(i,k,io1,0,lat)
          fk(i,ico) = fk(i,ico)-fs(i,k,ico,0,lat)
        enddo ! i=lon0,lon1

!       write(6,"(/,'comp: lat=',i3,' k=',i2)") lat,k
!       write(6,"('fk(:,io1)=',/,(6e12.4))") fk(lonbeg:lonend,io1)
!       write(6,"('fk(:,ico)=',/,(6e12.4))") fk(lonbeg:lonend,ico)
!
! Complete calculation of rhs in fk:
        do i=lonbeg,lonend
          fk(i,io1) = expz(k)*(o1nm_smooth(k,i,lat)*dtx2inv-fk(i,io1)+
     |      hdo1(k,i,lat))
          fk(i,ico) = expz(k)*(conm_smooth(k,i,lat)*dtx2inv-fk(i,ico)+
     |      hdco(k,i,lat))
        enddo ! i=lonbeg,lonend

!       write(6,"(/,'comp: lat=',i3,' k=',i2)") lat,k
!       write(6,"('fk(:,io1)=',/,(6e12.4))") fk(lonbeg:lonend,io1)
!       write(6,"('fk(:,ico)=',/,(6e12.4))") fk(lonbeg:lonend,ico)

!
! fk is ok up to this point.
! In earlier version, periodic points for fk were taken here.
! For now, ignore periodic points.


!
! Lower boundary:
        if (k==lev0) then
          do m=1,2                                                ! DO 16
            do kk=1,2
              do i=lon0,lon1
                qk(i,io1,m) = qk(i,io1,m)+pk(i,io1,kk)*b(i,kk,m) 
                qk(i,ico,m) = qk(i,ico,m)+pk(i,ico,kk)*b(i,kk,m) 
              enddo ! i=lon0,lon1
            enddo ! kk=1,2
          enddo ! m=1,2
          do m=1,2
            do i=lon0,lon1
              fk(i,io1) = fk(i,io1)-pk(i,io1,m)*fb(i,m)
              fk(i,ico) = fk(i,ico)-pk(i,ico,m)*fb(i,m)
              pk(i,:,m) = 0.
            enddo ! i=lon0,lon1
          enddo ! m=1,2

!         do m=1,2
!           write(6,"('comp lbc: m=',i2,' lat=',i2)") m,lat
!           write(6,"('qk(io1)=',/,(6e12.4))") qk(:,io1,m)
!           write(6,"('qk(ico)=',/,(6e12.4))") qk(:,ico,m)
!         enddo ! m=1,2
!         write(6,"('fk(io1)=',/,(6e12.4))") fk(:,io1)
!         write(6,"('fk(ico)=',/,(6e12.4))") fk(:,ico)

!
! Upper boundary:
        elseif (k==lev1-1) then
          do m=1,2
            do i=lon0,lon1
              qk(i,io1,m) = qk(i,io1,m)+(1.+.5*ep(i,m,kp)*dz)/
     |          (1.-.5*ep(i,m,kp)*dz)*rk(i,io1,m)
              qk(i,ico,m) = qk(i,ico,m)+(1.+.5*ep(i,m,kp)*dz)/
     |          (1.-.5*ep(i,m,kp)*dz)*rk(i,ico,m)
              rk(i,:,m) = 0.
            enddo ! i=lon0,lon1
          enddo ! m=1,2

!         do m=1,2
!           write(6,"('comp ubc: m=',i2,' lat=',i2)") m,lat
!           write(6,"('qk(io1)=',/,(6e12.4))") qk(:,io1,m)
!           write(6,"('qk(ico)=',/,(6e12.4))") qk(:,ico,m)
!         enddo ! m=1,2

        endif ! lbc or ubc
!
! QK=ALFAK=QK-PK*GAMA(K-1)
        do m=1,2                                          ! DO 18
          do kk=1,2
            do i=lon0,lon1

!             write(6,"('comp: i=',i2,' kk=',i2,' m=',i2,' k=',i2,
!    |          ' lat=',i2)") i,kk,m,k,lat
!             write(6,"('qk=',e12.4,' pk=',e12.4,' gama=',e12.4))") 
!    |          qk(i,io1,m),pk(i,io1,kk),gama(i,k,kk,m)

              qk(i,io1,m) = qk(i,io1,m)-pk(i,io1,kk)*gama(i,k,kk,m)
              qk(i,ico,m) = qk(i,ico,m)-pk(i,ico,kk)*gama(i,k,kk,m)
            enddo ! i=lon0,lon1
          enddo ! kk=1,2
        enddo ! m=1,2

! Testing will not work here until gama is incremented below..
!       do m=1,2
!         write(6,"('comp: m=',i2,' k=',i2,' lat=',i2)") m,k,lat
!         write(6,"('qk(io1)=',/,(6e12.4))") qk(:,io1,m)
!         write(6,"('qk(ico)=',/,(6e12.4))") qk(:,ico,m)
!       enddo ! m=1,2

!
! WKS1=DET(ALFA)
        do i=lon0,lon1
          wks1(i) = qk(i,1,1)*qk(i,2,2)-qk(i,1,2)*qk(i,2,1)
        enddo ! i=lon0,lon1
!
! WKM1=ALFAI
        do m=1,2
          do i=lon0,lon1
            wkm1(i,io1,m) = (delta(io1,m)*qk(i,ico,ico)-
     |        (1.-delta(io1,m))*qk(i,io1,m))/wks1(i)
            wkm1(i,ico,m) = (delta(ico,m)*qk(i,io1,io1)-
     |        (1.-delta(ico,m))*qk(i,ico,m))/wks1(i)
          enddo ! i=lon0,lon1
        enddo ! m=1,2
!
! WKV1=FK-PK*Z(K)
        do i=lon0,lon1
          wkv1(i,io1) = fk(i,io1)
          wkv1(i,ico) = fk(i,ico)
        enddo ! i=lon0,lon1
!
! GAMA(K+1)=ALFAI*RK
        do m=1,2
          do i=lon0,lon1
            gama(i,k+1,io1,m) = 0.
            gama(i,k+1,ico,m) = 0.
            wkv1(i,io1) = wkv1(i,io1)-pk(i,io1,m)*zz(i,k,m)
            wkv1(i,ico) = wkv1(i,ico)-pk(i,ico,m)*zz(i,k,m)
          enddo ! i=lon0,lon1
          do kk=1,2
            do i=lon0,lon1
              gama(i,k+1,io1,m) = gama(i,k+1,io1,m)+wkm1(i,io1,kk)*
     |          rk(i,kk,m)
              gama(i,k+1,ico,m) = gama(i,k+1,ico,m)+wkm1(i,ico,kk)*
     |          rk(i,kk,m)
            enddo ! i=lon0,lon1
          enddo ! kk=1,2
        enddo ! m=1,2
!
! Z(K+1)=WKM1*WKV1
        do i=lon0,lon1
          zz(i,k+1,:) = 0.
        enddo ! i=lon0,lon1
        do m=1,2
          do i=lon0,lon1
            zz(i,k+1,io1) = zz(i,k+1,io1)+wkm1(i,io1,m)*wkv1(i,m)
            zz(i,k+1,ico) = zz(i,k+1,ico)+wkm1(i,ico,m)*wkv1(i,m)
          enddo ! i=lon0,lon1
        enddo ! m=1,2

!       do m=1,2
!         write(6,"('comp: m=',i2,' k=',i2,' lat=',i2)") m,k,lat
!         write(6,"('gama(k+1,io1)=',/,(6e12.4))") gama(:,k+1,io1,m)
!         write(6,"('gama(k+1,ico)=',/,(6e12.4))") gama(:,k+1,ico,m)
!       enddo
!       write(6,"('zz(k+1,io1)=',/,(6e12.4))") zz(:,k+1,io1)
!       write(6,"('zz(k+1,ico)=',/,(6e12.4))") zz(:,k+1,ico)

!
! End main pressure loop:
      enddo levloop ! k=lev0,lev1-1
!
! Save diagnostics:
!     real,dimension(lev0:lev1,lon0:lon1,2) :: zz_ki
!     real,dimension(lev0:lev1,lon0:lon1,2,2) :: gama_ki
!
      do k=lev0,lev1
        zz_ki(k,:,io1) = zz(:,k,io1)
        zz_ki(k,:,ico) = zz(:,k,ico)
        do m=1,2
          gama_ki(k,:,io1,m) = gama(:,k,io1,m)
          gama_ki(k,:,ico,m) = gama(:,k,ico,m)
        enddo
      enddo
      call addfld('ZZ_O',' ',' ',zz_ki(:,:,io1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('ZZ_CO',' ',' ',zz_ki(:,:,ico),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('GAMAOM1',' ',' ',gama_ki(:,:,io1,1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('GAMAOM2',' ',' ',gama_ki(:,:,io1,2),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('GAMACOM1',' ',' ',gama_ki(:,:,ico,1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('GAMACOM2',' ',' ',gama_ki(:,:,ico,2),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Set upper boundary to zero:
      do i=lon0,lon1
        o1_upd(lev1,i,lat) = 0.
        co_upd(lev1,i,lat) = 0.
        upd(lev1,i,:) = 0.
      enddo ! i=lon0,lon1
!
! Downward sweep:
      do k=lev1-1,lev0,-1
        do isp=io1,ico
          do i=lon0,lon1
            upd(k,i,isp) = zz(i,k+1,isp)
          enddo ! i=lon0,lon1
          do m=1,2
            do i=lon0,lon1
              upd(k,i,isp) = upd(k,i,isp)-gama(i,k+1,isp,m)*
     |                       upd(k+1,i,m)
            enddo
          enddo ! m=1,2
        enddo ! isp=io1,ico
      enddo ! k=lev1-1,lev0,-1
!
! Transfer to output arrays:
      do k=lev0,lev1
        o1_upd(k,lon0:lon1,lat) = upd(k,:,io1)
        co_upd(k,lon0:lon1,lat) = upd(k,:,ico)
!  ** Overwriting update to test cooling/heating terms **
!       o1_upd(k,lon0:lon1,lat) = oref(k,lon0:lon1,lat)
!       co_upd(k,lon0:lon1,lat) = coref(k,lon0:lon1,lat)
      enddo ! k=lev0,lev1
!     do k=lev0,lev1
!  ** Overwriting update to test cooling/heating terms near the poles **
!       if (lat <= 14) then  !~56 degrees N/S
!       	o1_upd(k,lon0:lon1,lat) = o1_upd(k,lon0:lon1,lat)*1.05
!       elseif (lat >= 59) then  !~56 degrees N/S
!       	o1_upd(k,lon0:lon1,lat) = o1_upd(k,lon0:lon1,lat)*1.05
!       else
!       	o1_upd(k,lon0:lon1,lat) = o1_upd(k,lon0:lon1,lat)
!       endif
!     enddo ! k=lev0,lev1
!
! Upper boundary:
! kp is carried forward from the last iteration of levloop above.
      do i=lon0,lon1
        o1_upd(lev1,i,lat) = 
     |    (1.+.5*ep(i,io1,kp)*dz)/
     |    (1.-.5*ep(i,io1,kp)*dz)*o1_upd(lev1-1,i,lat)
        co_upd(lev1,i,lat) = 
     |    (1.+.5*ep(i,ico,kp)*dz)/
     |    (1.-.5*ep(i,ico,kp)*dz)*co_upd(lev1-1,i,lat)
      enddo ! i=lon0,lon1

      call addfld('O_SOLV',' ',' ',o1_upd(:,lon0:lon1,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('CO_SOLV',' ',' ',co_upd(:,lon0:lon1,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! End latitude scan:
      enddo ! lat=lat0,lat1
!
!
! Filter the new composition species:
!
! Set kut for wave filtering according to dlat (2.5 or 5.0):
!     call set_wave_filter(36,kut_5,nlat,kutt)
!
! Fourier smoothing of O and CO:
!     call filter_comp(o1_upd,lev0,lev1,lon0,lon1,lat0,lat1,kutt)
!     call filter_comp(co_upd,lev0,lev1,lon0,lon1,lat0,lat1,kutt)
      call filter_comp(o1_upd,lev0,lev1,lon0,lon1,lat0,lat1,kut)
      call filter_comp(co_upd,lev0,lev1,lon0,lon1,lat0,lat1,kut)
!
! Resume latitude scan:
      do lat=lat0,lat1

        call addfld('O_FILT',' ',' ',o1_upd(:,lon0:lon1,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('CO_FILT',' ',' ',co_upd(:,lon0:lon1,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Time smoothing:
        do i=lon0,lon1 
          do k=lev0,lev1
            o1nm_upd(k,i,lat) = dtsmooth_div2*(o1_nm(k,i,lat)+
     |        o1_upd(k,i,lat)) + dtsmooth*o1(k,i,lat)
            conm_upd(k,i,lat) = dtsmooth_div2*(co_nm(k,i,lat)+
     |        co_upd(k,i,lat)) + dtsmooth*co(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1 

!
      enddo ! lat=lat0,lat1

#ifdef MPI
!
! Periodic points:
!     call mp_periodic_f3d(o1_upd(:,lon0:lon1,lat0-1:lat1+1),
!    |  lev0,lev1,lon0,lon1,lat0-1,lat1+1)
!     call mp_periodic_f3d(co_upd(:,lon0:lon1,lat0-1:lat1+1),
!    |  lev0,lev1,lon0,lon1,lat0-1,lat1+1)
#endif
!
! Insure non-negative O,CO:
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
!           if (co_upd(k,i,lat) < small) co_upd(k,i,lat) = small
!           if (conm_upd(k,i,lat) < small) conm_upd(k,i,lat) = small
            if (co_upd(k,i,lat) < 1.0E-6) co_upd(k,i,lat) = 1.e-06
            if (conm_upd(k,i,lat) < 1.0E-6) conm_upd(k,i,lat) = 1.E-06
            if (o1_upd(k,i,lat) < small) o1_upd(k,i,lat) = small
            if (o1nm_upd(k,i,lat) < small) o1nm_upd(k,i,lat) = small
!           if (1.-small-o1_upd(k,i,lat)-co_upd(k,i,lat) < 0.) then
!             o1_upd(k,i,lat) = o1_upd(k,i,lat)*((1.-small)/
!    |                         (o1_upd(k,i,lat)+co_upd(k,i,lat)))
!             co_upd(k,i,lat) = co_upd(k,i,lat)*((1.-small)/
!    |                         (o1_upd(k,i,lat)+co_upd(k,i,lat)))
!           endif
!           if (1.-small-o1nm_upd(k,i,lat)-conm_upd(k,i,lat) < 0.) then
!             o1nm_upd(k,i,lat) = o1nm_upd(k,i,lat)*((1.-small)/
!    |                           (o1nm_upd(k,i,lat)+conm_upd(k,i,lat)))
!             conm_upd(k,i,lat) = conm_upd(k,i,lat)*((1.-small)/
!    |                           (o1nm_upd(k,i,lat)+conm_upd(k,i,lat)))
!           endif
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1

!       call addfld('ONM_OUT',' ',' ',o1nm_upd(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('CONM_OUT',' ',' ',conm_upd(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)

        call addfld('O_FINAL',' ',' ',o1_upd(:,lon0:lon1,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('CO_FINAL',' ',' ',co_upd(:,lon0:lon1,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('OMFINAL',' ',' ',o1nm_upd(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('COMFINAL',' ',' ',conm_upd(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
      enddo ! lat=lat0,lat1
      end subroutine comp
!-----------------------------------------------------------------------
      subroutine advecl(o1,co,un,vn,o1_advec,co_advec,
     |  lev0,lev1,lon0,lon1,lat0,lat1,lat)
!
! Horizontal advection for O,CO. 
! In previous versions, this was sub advecl (inline.F), called from 
!   k-loop in comp.F. Here it is called from beginning of comp.F at 
!   all latitudes and saved for later use (fk) inside comp.F k-loop.
! o,co,un,vn already have i-2,i-1,i+1,i+2, and j-1,j-2,j+1,j+2 for 
!   finite differencing.
!
      use cons_module,only: dlamda_2div3 ,dlamda_1div12, dphi_2div3,
     |  dphi_1div12,re_inv,racs
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  o1,co,un,vn
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(out) ::
     |  o1_advec, co_advec
!
! Local:
      integer :: k,i
!
      do k=lev0,lev1-1
        do i=lon0,lon1
          o1_advec(k,i,lat) = .5*racs(lat)*
     |      (dlamda_2div3*(o1(k,i+1,lat)-o1(k,i-1,lat))*
     |      (un(k,i+1,lat)+un(k,i-1,lat))-
     |      dlamda_1div12*(o1(k,i+2,lat)-o1(k,i-2,lat))*
     |      (un(k,i+2,lat)+un(k,i-2,lat)))+
     |      .5*re_inv*
     |      (dphi_2div3*(o1(k,i,lat+1)-o1(k,i,lat-1))*
     |      (vn(k,i,lat+1)+vn(k,i,lat-1))-
     |      dphi_1div12*(o1(k,i,lat+2)-o1(k,i,lat-2))*
     |      (vn(k,i,lat+2)+vn(k,i,lat-2)))

          co_advec(k,i,lat) = .5*racs(lat)*
     |      (dlamda_2div3*(co(k,i+1,lat)-co(k,i-1,lat))*
     |      (un(k,i+1,lat)+un(k,i-1,lat))-
     |      dlamda_1div12*(co(k,i+2,lat)-co(k,i-2,lat))*
     |      (un(k,i+2,lat)+un(k,i-2,lat)))+
     |      .5*re_inv*
     |      (dphi_2div3*(co(k,i,lat+1)-co(k,i,lat-1))*
     |      (vn(k,i,lat+1)+vn(k,i,lat-1))-
     |      dphi_1div12*(co(k,i,lat+2)-co(k,i,lat-2))*
     |      (vn(k,i,lat+2)+vn(k,i,lat-2)))
        enddo ! i=lon0,lon1
      enddo ! k=lev0,lev1-1

      end subroutine advecl
!-----------------------------------------------------------------------
      subroutine filter_comp(fout,lev0,lev1,lon0,lon1,lat0,lat1,kut)
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
#include <VT.inc>
#endif
!
! Local:
      integer :: i,j,nlevs,nlons,nlats
      real :: fik(nlonp4,lev0:lev1),fkij(lev0:lev1,nlonp4,lat0:lat1)
      real :: fmin,fmax
!
#ifdef VT
!     code = 131 ; state = 'filter_comp' ; activity='Filtering'
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
!     code = 131 ; state = 'filter_comp' ; activity='Filtering'
      call vtend(131,ier)
#endif
      end subroutine filter_comp
