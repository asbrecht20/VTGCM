!
      subroutine addiag(tn,o1,co,n2,vn,vc,barm,co2,z,
     |  lon0,lon1,lev0,lev1,lat0,lat1)
!
! Calculate needed terms vc, barm, and Z:
!  -- latest version: BF (8/1/07): add 75.0 km to lb heights
!  -- latest version: SWB (8/2/07): add 72.0 km to lb heights
!  -- latest version: SWB (11/27/07): add 70.0 km to lb heights
!  -- latest version: SWB (01/30/08): add updated gmlbc_oco to barm
!  -- latest version: SWB (01/30/08): add 68.0 km to lb heights
!  -- latest version: SWB (02/06/08): add 69.0 km to lb heights
!
      use cons_module,only: cs,dz,dzgrav,freq_semidi,dt,gmlbc_oco,
     |  gmlbc_oco_dbl,psin2,rmassinv_o1,rmassinv_co,rmassinv_co2,
     |  boltz,expz,p0
      use lbc,only: z_lbc
      use init_module,only: iter
      use fields_module,only: tlbc
      use addfld_module,only: addfld
      use input_module,only: oxvgcm_ncfile
      implicit none
!
! Input args:
      integer,intent(in) :: lon0,lon1,lev0,lev1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn, ! neutral temperature (deg K)
     |  o1, ! atomic oxygen (mmr)
     |  co, ! carbon monoxide (mmr)
     |  n2, ! nitrogen dioxide (mmr)
     |  vn  ! meridional wind velocity (cm/s)
!
! Output args:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  vc  ,
     |  barm,
     |  co2 ,
     |  z
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
! Local:
      integer :: k,i,j,nlevs,nlons,nlats,ier
      real ::
     |  barm1(lon0:lon1),
     |  w1   (lev0:lev1,lon0:lon1),
     |  xnmbar   (lev0:lev1,lon0:lon1),
     |  xtot     (lev0:lev1,lon0:lon1,lat0:lat1)
      complex :: expt
      real :: fmin,fmax,pso,psco
!
#ifdef VT
!     code = 116 ; state = 'addiag' ; activity='ModelCode'
      call vtbegin(116,ier)
#endif
      nlevs = lev1-lev0+1
      nlons = lon1-lon0+1
      nlats = lat1-lat0+1
!
! Latitude scan:
      do j=lat0,lat1

        call addfld('TN_diag',' ',' ',tn(:,lon0:lon1,j),
     |    'lev',lev0,lev1,'lon',lon0,lon1,j)
        call addfld('VN_diag',' ',' ',vn(:,lon0:lon1,j),
     |    'lev',lev0,lev1,'lon',lon0,lon1,j)
        call addfld('O1_diag',' ',' ',o1(:,lon0:lon1,j),
     |    'lev',lev0,lev1,'lon',lon0,lon1,j)
        call addfld('CO_diag',' ',' ',co(:,lon0:lon1,j),
     |    'lev',lev0,lev1,'lon',lon0,lon1,j)
        call addfld('N2_diag',' ',' ',n2(:,lon0:lon1,j),
     |    'lev',lev0,lev1,'lon',lon0,lon1,j)
!
! vc = cos(phi)*v
        do i=lon0,lon1
          do k=lev0,lev1
            vc(k,i,j) = cs(j)*vn(k,i,j)
          enddo
        enddo
        call addfld('VC_diag',' ',' ',vc(lev0:lev1-1,lon0:lon1,j),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,j)
!
! barm = mean molecular weight (k+1/2):
!
        do i=lon0,lon1
          do k=lev0,lev1

! vtgcm: PSCO2=(1.-F(I,NPS1K)-F(I,NPS2K)-PSIN2(I,K,J+2))
! vtgcm co2 = 1.-o-co-n2:
            co2(k,i,j) = (1.-o1(k,i,j)-co(k,i,j)-n2(k,i,j))
!
! co2 < 0. only at lev1:
            if (co2(k,i,j) < 0.) then
!             write(6,"('>>> addiag warning: co2(k,i,j) < 0: ',
!    |          ' k,i,j=',3i4,' co2=',e12.4)") k,i,j,co2(k,i,j)
              co2(k,i,j) = 0.
            endif

!! vtgcm:
!!     F(I,NMSK)=1./(F(I,NPS1K)/RMASS(1)+(F(I,NPS2K)+PSIN2(I,K,J+2))/
!!    & RMASS(2)+PSCO2/RMASS(3))
!! (NPS1=o, NPS2=co, rmass = /o,co,co2/)

!
! barm = 1/(o/16+(co+n2)/28)+co2/44
            barm(k,i,j) = 1./ 
     |        (o1(k,i,j)*rmassinv_o1+(co(k,i,j)+n2(k,i,j))*
     |         rmassinv_co+co2(k,i,j)*rmassinv_co2)
          enddo
        enddo
        call addfld('BARM1',' ',' ',barm(:,lon0:lon1,j),
     |    'lev',lev0,lev1,'lon',lon0,lon1,j)
!
! xnmbari = p0*e(-z)*barm/kT at interfaces:
        do i=lon0,lon1
          do k=lev0,lev1
            xnmbar(k,i) = p0*expz(k)*barm(k,i,j)/
     |        (boltz*tn(k,i,j))
          enddo
        enddo
        do i=lon0,lon1
          do k=lev0,lev1
!           xtot(k,i,j) = xnmbar(k,i)/barm(k,i,j)
            xtot(k,i,j) = p0*expz(k)/(boltz*tn(k,i,j))
          enddo
        enddo
        call addfld('XTOT',' ',' ',xtot(:,lon0:lon1,j),
     |    'lev',lev0,lev1,'lon',lon0,lon1,j)

!
! 9/19/06 btf: Currently, tgcmproc_f90 makes CO2 a derived field,
!   dependent on O, CO, and N2, so do not use "CO2" name here:
!
        call addfld('CO2_diag',' ',' ',co2(:,lon0:lon1,j),
     |    'lev',lev0,lev1,'lon',lon0,lon1,j)
!
! barm1 = barm(k=0) (linear extrapolation)
!
        do i=lon0,lon1
          barm1(i) = 1.5*barm(lev0,i,j)-0.5*barm(lev0+1,i,j)
        enddo
!
! barm(k) = 0.5*(barm(k+1/2)+barm(k-1/2)), k = kmaxp1,2,1
!
        do i=lon0,lon1
          do k=lev1,lev0+1,-1
            barm(k,i,j) = 0.5*(barm(k,i,j)+barm(k-1,i,j))
          enddo
        enddo
!
! barm(1) = barm1
!
        do i=lon0,lon1
          barm(lev0,i,j) = barm1(i)
        enddo
        call addfld('BARM2',' ',' ',barm(:,lon0:lon1,j),
     |    'lev',lev0,lev1,'lon',lon0,lon1,j)
!
! Lower boundary of barm (from mtgcm16):
!    |  gmlbc_oco    = 2.*0,                   ! vtgcm only       C(46)
!    |  gmlbc_oco_dbl= 2.*2.64E-04,            ! vtgcm only       C(47)
!    |  psin2        = 0.0265*(28./44.),       ! vtgcm only       GMEAN
! 
        do i=lon0,lon1
          pso = gmlbc_oco
          psco = gmlbc_oco_dbl
          co2(lev0,i,j) = (1.-pso-psco-psin2)
          barm(lev0,i,j) = 1./(pso*rmassinv_o1+(psco+psin2)*rmassinv_co+
     |      co2(lev0,i,j)*rmassinv_co2)
        enddo ! i=lon0,lon1
!
! Calculate z:
!
! w1 = barm
        do i=lon0,lon1
          do k=lev0,lev1-1
            w1(k,i) = (barm(k,i,j)+barm(k+1,i,j))*0.5
          enddo
        enddo
!
! w1 = tn/w1 (old model comment: s1=s2/s1=(t+t0)/m)
        do i=lon0,lon1
          do k=lev0,lev1-1
            w1(k,i) = tn(k,i,j)/w1(k,i)
          enddo
        enddo
        call addfld('W1a_diag',' ',' ',w1(lev0:lev1-1,:),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,j)
!
! w1=(ds*r/g)*w1
        do i=lon0,lon1
          do k=lev0,lev1-1
            w1(k,i) = (dz/dzgrav) * w1(k,i)
          enddo
        enddo
        call addfld('W1b_diag',' ',' ',w1(lev0:lev1-1,:),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,j)
!
! Lower boundary of geopotential from lbc.F:
        do i=lon0,lon1
          z(lev0,i,j) = z_lbc(i,j)
        enddo
!       call addfld('Z_LBC',' ',' ',z_lbc(lon0:lon1,j),
!    |    'lon',lon0,lon1,'lat',j,j,0)
!
! z(k+1)=s1(k)+z(k)
!
        do i=lon0,lon1
          do k=lev0,lev1-1
            z(k+1,i,j) = w1(k,i)+z(k,i,j)
          enddo
        enddo
!
! 9/06 btf: add 94 km to Z:
!       z(:,:,j) = z(:,:,j)+9.4e6 ! cm
! 7/30/07 btf: add 75 km to Z (new lbc zp -16)
!       z(:,:,j) = z(:,:,j)+7.5e6 ! cm
! 8/2/07 swb: add 72 km to Z (new lbc zp -16)
!       z(:,:,j) = z(:,:,j)+7.2e6 ! cm
! 11/27/07 swb: add 70 km to Z (new lbc zp -16)
!       z(:,:,j) = z(:,:,j)+7.0e6 ! cm
! 01/30/08 swb: add 68.0 km to Z (new lbc zp -16)
! (only if not nudging)
        if (len_trim(oxvgcm_ncfile) == 0)
     |    z(:,:,j) = z(:,:,j)+6.9e6 ! cm
!
        call addfld('Z_diag','Geopotential from addiag',
     |    'cm',z(:,lon0:lon1,j),'ilev',lev0,lev1,'lon',lon0,lon1,j)
!
!       do k=lev0,lev1-1
!         write(6,"('addiag: k=',i3,' j=',i3,' z(k,lon0:lon1,j)=',
!    |      /,(6e12.4))") k,j,z(k,lon0:lon1,j)
!       enddo ! k=lev0,lev1-1
!
! End latitude scan:
      enddo ! j=lat0,lat1
!
#ifdef VT
!     code = 116 ; state = 'addiag' ; activity='ModelCode'
      call vtend(116,ier)
#endif
      end subroutine addiag
