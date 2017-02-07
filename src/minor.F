!
      subroutine minor(tn,difk,o1,co,co2,n2,fcomp,fcomp_tm1,fcomp_out,
     |  fcomp_tm1_out,sloss,sprod,flbc,fubc,rmx,phix,alfax,
     |  lev0,lev1,lon0,lon1,lat0,lat1,idebug)
!
! Advance minor species fcomp by one time step.
! This is called from comp_o2 and comp_n4s.
! Last modified by SWB:  2/15/07  (vtgcm changes)
! Last modified by SWB:  5/06/08  call trsolv changes (for lev1-1)
! Last modified by SWB:  5/14/08  psin2 = 0.0277
! Last modified by SWB:  5/15/08  small = 1.0E-14 
! Last modified by SWB:  5/15/08  small = 1.0E-18 
! Last modified by SWB:  5/19/08  small = 1.0E-24 
! Last modified by SWB:  5/30/08  psin2 = 0.064
!
      use params_module,only: nlat,nlon,nlonp4,dlat
      use lbc,only: b,fb ! b(nlonp4,2,2), fb(nlonp4,2)
!     vtgcm changes: SWB
      use cons_module,only: rmassinv_o2,rmassinv_o1,rmassinv_n2,
     |  rmassinv_co2,rmassinv_co,rmass_co2,rmass_co,dzp,
     |  boltz,p0,expz,rmass_o2,rmass_o1,rmass_n2,hor,dtr,pi,shapiro,
     |  expzmid_inv,expzmid,dtx2inv,grav,avo,
     |  kut,dtsmooth,dtsmooth_div2
      use init_module,only: glat
      use input_module,only: difhor
      use fields_module,only: w,itc,itp
#ifdef MPI
      use mpi_module,only: mp_bndlons_f3d, mp_periodic_f3d
#endif
!
! VT means vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
      implicit none
!
! Args:
      integer :: lev0,lev1,lon0,lon1,lat0,lat1,idebug
!
! Input fields are at full task subdomain:
!     vtgcm changes: SWB
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn,       ! neutral temperature (deg K)
     |  difk,     ! eddy diffusion (eddy_flds)
     |  o1,       ! atomic oxygen (mmr)
     |  co,       ! carbon monoxide (mmr)
     |  co2,      ! carbon dioxide (mmr)
     |  n2,       ! molecular nitrogen (mmr)
     |  fcomp,    ! input species (for finite diffs in advec)
     |  fcomp_tm1 ! input species at time n-1 
!
! flbc(:,1)=A, flbc(:,2)=B, flbc(:,3)=C define lower boundary condition, 
!   where: A*DPSX/DZ + B*PSX + C = 0.
! Boundary conditions are allocated at task subdomains, without ghost cells.
!
      real,dimension(lon0:lon1,3,lat0:lat1),intent(in) :: flbc
!
! fubc = diffusive upward number flux at upper boundary.
      real,dimension(lon0:lon1,lat0:lat1),intent(in) :: fubc
!
      real,intent(in) :: 
     |  phix(3), ! diffusion vector
     |  rmx,     ! molecular weight of fcomp minor species
     |  alfax    ! thermal diffusion coefficient
!
! Input production and loss:
! sloss: sx/n(x) where sx is portion of number density source
!   proportional to n(x), the minor species number density.
! sprod: s0, portion of number density source independent of n(x)
!
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(in) ::
     |  sloss, sprod
!
! Output minor species (allocated at full subdomains):
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  fcomp_out,     ! output species 
     |  fcomp_tm1_out  ! output species at time n-1
!
! Local:
      integer :: k,i,lonbeg,lonend,lat,ier
!
      integer :: nk,nkm1,nlevs
!     vtgcm changes: SWB
!     real,parameter :: small=1.e-12, tau=9.25e+3, t00=273.,
!     |  psin2b = 0.0261
!     real,parameter :: small=1.e-18, tau=9.25e+3, t00=273.,
!    |  psin2b = 0.044
      real,parameter :: small=1.e-24, tau=9.25e+3, t00=273.,
     |  psin2b = 0.064
      real :: salfa12,salfa21,salfax1,salfax2,rlat,dfac
!     vtgcm changes: SWB
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  hadvec,   ! horizontal advection (output of sub advec) (s13)
!    |  do2dz,    ! do2/dz(k)    (s6)
!    |  do1dz,    ! do1/dz(k)    (s7)
     |  do1dz,    ! do1/dz(k)    (s6)
     |  dcodz,    ! dco/dz(k)    (s7)
!    |  pso2,     ! o2           (s8)
!    |  pso1,     ! o1           (s9)
     |  pso1,     ! o1 mmr       (s8)
     |  psco,     ! co mmr       (s9)
     |  dmdz,     ! dm/dz(k)     (s10)
     |  xmbar_k,  ! mbar (k)     (s11)
     |  xmbar_kh, ! mbar (k+1/2) (s12)
!    |  xn2,      ! N2 (mmr)
     |  psn2,     ! N2 (mmr)
     |  psco2,    ! co2 (mmr)
     |  tni,      ! tn (k)       (s5)
     |  s0prod,   ! sprod*mx/nmbar (s15)
     |  alfa11,alfa12,alfa21,alfa22, ! (s1,s2,s3,s4)
     |  ex,ax,    ! (s12,s11)
     |  thdiff,   ! thermal diffusion term (s12)
     |  p_coef,q_coef,r_coef,f_rhs ! coefficients for tridiagonal solver
!
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  ftm1_smooth  ! time-smoothed field at time n-1 (s10)
!
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0:lat1) ::
     |  ftm1_jdif ! lat diffs (includes lon boundaries) (s8)
!
      integer :: jm2,jm1,j0,jp1,jp2 ! lat-2,lat-1,lat,lat+1,lat+2
!
!     vtgcm changes: SWB
      real,dimension(lon0:lon1) :: 
     |  xmbari,   ! (t6)
!    |  bo2,bo1,  ! (t4,t5)
     |  bo1,bco,  ! (t4,t5)
     |  dfactor   ! (t7) (formerly output of sub dfact)
!
      real :: phi(2,3)
!
#ifdef VT
!     code = 114 ; state = 'minor' ; activity='ModelCode'
      call vtbegin(114,ier)
#endif
!
      if (idebug > 0) write(6,"('Enter minor')")

!     vtgcm changes: SWB
      phi(1,:) = (/0.   , 1.31,  1.026 /)
      phi(2,:) = (/2.296, 0.  ,  1.600/)

      salfa12=phi(1,2)-phi(1,3)
      salfa21=phi(2,1)-phi(2,3)
      salfax1=phix(1)-phix(3)
      salfax2=phix(2)-phix(3)

      nk = lev1-lev0+1
      nkm1 = nk-1
      nlevs = nk

      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = lon1-2 ! nlonp4==nlonp4
!
! Latitude scan over task subdomain:
      do lat=lat0,lat1
!       if (idebug > 0) 
!    |    write(6,"('Minor first latitude scan: lat=',i3)") lat
        jm2 = lat-2
        jm1 = lat-1
        j0  = lat
        jp1 = lat+1
        jp2 = lat+2
!
! Shapiro time smoother:
!
! ftm1_jdif = 4th order diffs in latitude at time n-1:
! f = f(j)-shapiro * ( f(j+2)+f(j-2) - 4.*(f(j+1)+f(j-1)) + 6.*f(j) )
!
        do i=lon0,lon1
          do k=lev0,lev1
            ftm1_jdif(k,i,lat) = fcomp_tm1(k,i,j0) - shapiro * 
     |           (fcomp_tm1(k,i,jp2) + fcomp_tm1(k,i,jm2)  - 
     |        4.*(fcomp_tm1(k,i,jp1) + fcomp_tm1(k,i,jm1)) + 
     |        6.*fcomp_tm1(k,i,j0))
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1

      enddo ! lat=lat0,lat1

!     if (idebug > 0) 
!    |  write(6,"('Minor after first lat scan: call bndlons')")

!------------------------- End first latitude scan ---------------------
#ifdef MPI
!
! Boundary longitudes for ftm1_jdif:
! 10/24/02 bf: changed following call:
!     call mp_bndlons_f3d(ftm1_jdif,nlevs,lon0,lon1,lat0-2,lat1+2,1)
      call mp_bndlons_f3d(ftm1_jdif,nlevs,lon0,lon1,lat0,lat1,1)

!     if (idebug > 0) 
!    |  write(6,"('Minor after bndlons_f3d')")
#endif
!----------------------- Begin second latitude scan --------------------
      do lat=lat0,lat1
!       if (idebug > 0) 
!    |    write(6,"('Minor second latitude scan: lat=',i3)") lat
!
! ftm1_smooth = zonally-smoothed field at time n-1:
! f = f(i)-shapiro * ( f(i+2)+f(i-2) - 4.*(f(i+1)+f(i-1)) + 6.*f(i) )
      do i=lonbeg,lonend
        do k=lev0,lev1-1
          ftm1_smooth(k,i,lat) = ftm1_jdif(k,i,lat) - shapiro *
     |         (ftm1_jdif(k,i+2,lat) + ftm1_jdif(k,i-2,lat)  -
     |      4.*(ftm1_jdif(k,i+1,lat) + ftm1_jdif(k,i-1,lat)) +
     |      6.* ftm1_jdif(k,i,lat))
        enddo ! k=lev0,lev1-1
      enddo ! i=lonbeg,lonend
!
! Set periodic points to zero:
      if (lon0==1) ftm1_smooth(:,lon0:lon0+1,lat) = 0.
      if (lon1==nlonp4) ftm1_smooth(:,lon1-1:lon1,lat) = 0.

!     if (idebug > 0) write(6,"('Minor second lat scan after ',
!    |  ' ftm1_smooth: lat=',i3)") lat
!
! Horizontal advection (pass k vs i slices at full task subdomain 
! longitudes, and the 5 latitudes centered over the current latitude).
!
      call advec(fcomp(:,:,lat-2:lat+2),hadvec,lev0,lev1,lon0,lon1,
     |  lat)

!     if (idebug > 0) write(6,"('Minor second lat scan after ',
!    |  ' advec: lat=',i3)") lat
!
! Periodic points for advection:
!     call periodic_f2d(hadvec,lon0,lon1,nk)
!
! Set periodic points to zero:
      if (lon0==1) hadvec(:,lon0:lon0+1) = 0.
      if (lon1==nlonp4) hadvec(:,lon1-1:lon1) = 0.

!     if (idebug > 0) then
!       call addfsech('HADVEC',' ',' ',hadvec,lon0,lon1,nk,nkm1,lat)
!       call addfsech('SLOSS' ,' ',' ',sloss(:,lon0:lon1,lat),
!    |    lon0,lon1,nk,nkm1,lat)
!       call addfsech('SPROD' ,' ',' ',sprod(:,lon0:lon1,lat),
!    |    lon0,lon1,nk,nkm1,lat)
!     endif
!     vtgcm changes: SWB
      do i=lon0,lon1
!       Earth TIEGCM  -------------------------------------------------------
!       bo2(i) = b(i,1,1)*o2(lev0,i,lat)+b(i,1,2)*o1(lev0,i,lat)+fb(i,1) ! t4
!       bo1(i) = b(i,2,1)*o2(lev0,i,lat)+b(i,2,2)*o1(lev0,i,lat)+fb(i,2) ! t5
!       xmbari(i) = 1./(bo2(i)*rmassinv_o2+bo1(i)*rmassinv_o1+   ! t6
!    |    (1.-bo2(i)-bo1(i))*rmassinv_n2)
!       Venus TGCM ----------------------------------------------------------
        bo1(i) = b(i,1,1)*o1(lev0,i,lat)+b(i,1,2)*co(lev0,i,lat)+fb(i,1) ! t4
        bco(i) = b(i,2,1)*o1(lev0,i,lat)+b(i,2,2)*co(lev0,i,lat)+fb(i,2) ! t5
        xmbari(i) = 1./(bo1(i)*rmassinv_o1+bco(i)*rmassinv_co +
     |    psin2b*rmassinv_n2+ (1.-bo1(i)-bco(i)-psin2b)*rmassinv_co2)   ! t6
      enddo ! i=lon0,lon1

!     if (idebug > 0) write(6,"('Minor second lat scan after ',
!    |  ' xmbari: lat=',i3)") lat
!
! xmbar_kh = mbar(k+1/2):
!     vtgcm changes: SWB
      do i=lon0,lon1
        do k=lev0,lev1
!         xmbar_kh(k,i) = 1./(o2(k,i,lat)*rmassinv_o2+o1(k,i,lat)*
!    |      rmassinv_o1+xn2(k,i)*rmassinv_n2)
          xmbar_kh(k,i) = 1./(o1(k,i,lat)*rmassinv_o1+co(k,i,lat)*
     |      rmassinv_co+n2(k,i,lat)*rmassinv_n2+co2(k,i,lat)*
     |      rmassinv_co2)
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1

!     if (idebug > 0) write(6,"('Minor second lat scan after ',
!    |  ' xmbar_kh: lat=',i3)") lat
!
! Lower boundary:
!     vtgcm changes: SWB
      do i=lon0,lon1
        xmbar_k(lev0,i) = .5*(xmbari(i)+xmbar_kh(lev0,i))
        dmdz(lev0,i) = (xmbar_kh(lev0,i)-xmbari(i))/dzp
        pso1(lev0,i) = .5*(bo1(i)+o1(lev0,i,lat))  
        psco(lev0,i) = .5*(bco(i)+co(lev0,i,lat))  
        do1dz(lev0,i) = (o1(lev0,i,lat)-bo1(i))/dzp
        dcodz(lev0,i) = (co(lev0,i,lat)-bco(i))/dzp
      enddo ! i=lon0,lon1

!     if (idebug > 0) write(6,"('Minor second lat scan after ',
!    |  ' lbc: lat=',i3)") lat
!
! Levels 2 -> lev1:
!     vtgcm changes: SWB
      do i=lon0,lon1
        do k=lev0+1,lev1
          xmbar_k(k,i) = .5*(xmbar_kh(k,i)+xmbar_kh(k-1,i))
          dmdz(k,i) = (xmbar_kh(k,i)-xmbar_kh(k-1,i))/dzp
          pso1(k,i) = .5*(o1(k,i,lat)+o1(k-1,i,lat))
          psco(k,i) = .5*(co(k,i,lat)+co(k-1,i,lat))
          do1dz(k,i) = (o1(k,i,lat)-o1(k-1,i,lat))/dzp
          dcodz(k,i) = (co(k,i,lat)-co(k-1,i,lat))/dzp
        enddo ! k=lev0+1,lev1
      enddo ! i=lon0,lon1

!     if (idebug > 0) write(6,"('Minor second lat scan after ',
!    |  ' levels 2-lev1: lat=',i3)") lat
!
! tni = tn at interfaces:
      do i=lon0,lon1
        tni(lev0,i) = tn(lev1,i,lat) ! assumes bottom boundary of tn in top slot
        tni(lev1,i) = tn(lev1-1,i,lat)
        do k=lev0+1,lev1-1
          tni(k,i) = .5*(tn(k,i,lat)+tn(k-1,i,lat))
        enddo ! k=lev0+1,lev1-1
      enddo ! i=lon0,lon1

!     if (idebug > 0) write(6,"('Minor second lat scan after ',
!    |  ' tni: lat=',i3)") lat

!     if (idebug > 0) then
!       call addfsech('XMBAR_KH',' ',' ',xmbar_kh,lon0,lon1,nk,nkm1,lat) ! s12
!       call addfsech('XMBAR_K' ,' ',' ',xmbar_k ,lon0,lon1,nk,nkm1,lat) ! s11
!       call addfsech('DMDZ0'   ,' ',' ',dmdz    ,lon0,lon1,nk,nkm1,lat) ! s10
!       call addfsech('PSO1'    ,' ',' ',pso1    ,lon0,lon1,nk,nkm1,lat) ! s9
!       call addfsech('PSO2'    ,' ',' ',pso2    ,lon0,lon1,nk,nkm1,lat) ! s8
!       call addfsech('DO1DZ'   ,' ',' ',do1dz   ,lon0,lon1,nk,nkm1,lat) ! s7
!       call addfsech('DO2DZ'   ,' ',' ',do2dz   ,lon0,lon1,nk,nkm1,lat) ! s6
!       call addfsech('TNI'     ,' ',' ',tni     ,lon0,lon1,nk,nkm1,lat) ! s5
!     endif
 
      do i=lon0,lon1
        do k=lev0,lev1-1
          s0prod(k,i) = sprod(k,i,lat)*rmx*boltz*tn(k,i,lat)/  ! s15
     |      (p0*expz(k)*xmbar_kh(k,i))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     if (idebug > 0) write(6,"('Minor second lat scan after ',
!    |  ' s0prod: lat=',i3,' lat0,1=',2i4)") lat,lat0,lat1

!     vtgcm changes: SWB
      do i=lon0,lon1
        do k=lev0,lev1
!    -----------------------------------------------------
          alfa11(k,i) = -(phi(1,3)+salfa12*psco(k,i)) ! s1
          alfa12(k,i) = salfa12*pso1(k,i)             ! s2
          alfa21(k,i) = salfa21*psco(k,i)             ! s3
          alfa22(k,i) = -(phi(2,3)+salfa21*pso1(k,i)) ! s4
!    -----------------------------------------------------
!
!        if (idebug > 0)
!    |    write(6,"('Minor calc ex: k,i,lat=',3i4,' lat0,1=',2i4,
!    |      ' xmbar_k=',e12.4,' salfax1,2=',2e12.4,' alfa11,12,21,22=',
!    |      4e12.4)") k,i,lat,lat0,lat1,xmbar_k(k,i),salfax1,salfax2,
!    |      alfa11(k,i),alfa12(k,i),alfa21(k,i),alfa22(k,i)

          ex(k,i) =                                   ! s12
     |      ((salfax1*alfa22(k,i)-salfax2*alfa21(k,i))*(do1dz(k,i)-
     |      (1.-(rmass_o1+dmdz(k,i))/xmbar_k(k,i))*pso1(k,i))+
     |      (salfax2*alfa11(k,i)-salfax1*alfa12(k,i))*
     |      (dcodz(k,i)-(1.-(rmass_co+dmdz(k,i))/xmbar_k(k,i))*
     |      psco(k,i)))/(alfa11(k,i)*alfa22(k,i)-alfa12(k,i)*
     |      alfa21(k,i))+1.-(rmx+dmdz(k,i))/xmbar_k(k,i)
 
          dmdz(k,i) = dmdz(k,i)/xmbar_k(k,i)          ! s10
!
!         if (idebug > 0)
!    |    write(6,"('Minor calc ax: k,i,lat=',3i4,' lat0,1=',2i4,
!    |      ' phix(3)=',e12.4,' salfax1,2=',2e12.4,' pso1=',e12.4,
!    |      ' psco=',e12.4)") k,i,lat,lat0,lat1,phix(3),salfax1,salfax2,
!    |      pso1(k,i),psco(k,i)

          ax(k,i) = 
     |      -xmbar_k(k,i)/(tau*rmass_co2)*(t00/tni(k,i))**0.25/
     |      (phix(3)+salfax1*pso1(k,i)+salfax2*psco(k,i))

        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1

!     if (idebug > 0) write(6,"('Minor second lat scan after ',
!    |  ' ex,ax: lat=',i3)") lat

!     if (idebug > 0) then
!       call addfsech('S0PROD',' ',' ',s0prod,lon0,lon1,nk,nkm1,lat) ! s15
!       call addfsech('ALFA11',' ',' ',alfa11,lon0,lon1,nk,nk,lat) ! s1
!       call addfsech('ALFA12',' ',' ',alfa12,lon0,lon1,nk,nk,lat) ! s2
!       call addfsech('ALFA21',' ',' ',alfa21,lon0,lon1,nk,nk,lat) ! s3
!       call addfsech('ALFA22',' ',' ',alfa22,lon0,lon1,nk,nk,lat) ! s4
!       call addfsech('EX'    ,' ',' ',ex    ,lon0,lon1,nk,nk,lat) ! s12
!       call addfsech('DMDZ1' ,' ',' ',dmdz  ,lon0,lon1,nk,nk,lat) ! s10
!       call addfsech('AX0'   ,' ',' ',ax    ,lon0,lon1,nk,nk,lat) ! s11
!     endif
!
! thdiff = EX-ALFAX*D/DS(LN(T(TOT))  (Thermal diffusion term) (s12)
      do i=lon0,lon1
        do k=lev0+1,lev1-1
          thdiff(k,i) = ex(k,i)-alfax*(tni(k+1,i)-tni(k-1,i))/
     |      (2.*dzp*tni(k,i))
        enddo ! k=lev0+1,lev1-1
        thdiff(lev0,i) = ex(lev0,i)-alfax*(tni(lev0+1,i)-tni(lev0,i))/
     |    (dzp*tni(lev0,i))
        thdiff(lev1,i) = ex(lev1,i)-alfax*(tni(lev1,i)-tni(lev1-1,i))/
     |    (dzp*tni(lev1,i))
      enddo ! i=lon0,lon1

!     if (idebug > 0) write(6,"('Minor second lat scan after ',
!    |  ' thdiff: lat=',i3)") lat
!     if (idebug > 0) 
!    |  call addfsech('THDIFF0',' ',' ',thdiff,lon0,lon1,nk,nk,lat) ! s12
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

!     if (idebug > 0) write(6,"('Minor second lat scan after ',
!    |  ' dfactor: lat=',i3)") lat
!
! Set up coefficients for tridiagonal solver:
! (W (f4s(w)) is use-associated from fields module)
! p_coef (s1), q_coef (s2), r_coef (s3), and f_rhs (s4)
! (difk = eddy diffusion)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
! s1
          p_coef(k,i) = ax(k,i)/dzp*(1./dzp+.5*thdiff(k,i))-expz(k)*
     |      (expzmid_inv*difk(k,i,lat)*dfactor(i)*(1./dzp-.5*dmdz(k,i))+
     |      0.25*(w(k,i,lat,itc)+w(k+1,i,lat,itc)))/
     |      dzp
! s3
          r_coef(k,i) = ax(k+1,i)/dzp*(1./dzp-.5*thdiff(k+1,i))-expz(k)*
     |      (expzmid*difk(k+1,i,lat)*dfactor(i)*(1./dzp+.5*dmdz(k+1,i))-
     |      0.25*(w(k,i,lat,itc)+w(k+1,i,lat,itc)))/
     |      dzp
! s2
          q_coef(k,i) = 
     |      -(ax(k  ,i)/dzp*(1./dzp-.5*thdiff(k  ,i)) +
     |        ax(k+1,i)/dzp*(1./dzp+.5*thdiff(k+1,i)))+
     |      expz(k)*((expzmid_inv*difk(k,i,lat)*(1./dzp+.5*dmdz(k,i))+
     |      expzmid*difk(k+1,i,lat)*(1./dzp-.5*dmdz(k+1,i)))*
     |      dfactor(i)/dzp-sloss(k,i,lat)+dtx2inv)
! s4
          f_rhs(k,i) = expz(k)*(ftm1_smooth(k,i,lat)*
     |      dtx2inv-hadvec(k,i)+s0prod(k,i))
        enddo ! k=lev0,lev1-1
!
! Lower boundary (use flbc(:,1-3):
        q_coef(lev0,i) = q_coef(lev0,i)+p_coef(lev0,i)*(flbc(i,1,lat)+
     |    .5*flbc(i,2,lat)*dzp)/(flbc(i,1,lat)-.5*flbc(i,2,lat)*dzp)
        f_rhs(lev0,i) = f_rhs(lev0,i)-p_coef(lev0,i)*flbc(i,3,lat)*dzp/
     |    (flbc(i,1,lat)-.5*flbc(i,2,lat)*dzp)
        p_coef(lev0,i) = 0.
!
! Upper boundary (use fubc(:):
        p_coef(lev1,i) = 1.+.5*dzp*thdiff(lev1,i)
        q_coef(lev1,i) = p_coef(lev1,i)-2.
        r_coef(lev1,i) = 0.
        f_rhs (lev1,i) = -grav*rmx*fubc(i,lat)*dzp/(p0*ax(lev1,i)*avo)
      enddo ! i=lon0,lon1

!     if (idebug > 0) write(6,"('Minor second lat scan after ',
!    |  ' coeffs: lat=',i3)") lat
!     if (idebug > 0) then
!       call addfsech('AX1'    ,' ',' ',ax         ,lon0,lon1,nk,nk,lat) ! s11
!       call addfsech('DMDZ2'  ,' ',' ',dmdz       ,lon0,lon1,nk,nk,lat) ! s9
!       call addfsech('THDIFF1',' ',' ',thdiff     ,lon0,lon1,nk,nk,lat) ! s12
!       call addfsech('FSMOOTH',' ',' ',ftm1_smooth(:,lon0:lon1,lat),
!    |    lon0,lon1,nk,nkm1,lat) ! s10
!       call addfsech('W_OMEGA',' ',' ',w(lev0:lev1,lon0:lon1,
!    |    lat,itc),lon0,lon1,nk,nkm1,lat) ! s10
!
!       call addfsech('P_COEF' ,' ',' ',p_coef     ,lon0,lon1,nk,nk,lat) ! s1
!       call addfsech('R_COEF' ,' ',' ',r_coef     ,lon0,lon1,nk,nk,lat) ! s3
!       call addfsech('Q_COEF' ,' ',' ',q_coef     ,lon0,lon1,nk,nk,lat) ! s2
!       call addfsech('F_RHS'  ,' ',' ',f_rhs      ,lon0,lon1,nk,nk,lat) ! s4
!     endif
!
! Solve tridiagonal system (fcomp_out is allocated at full subdomain):
!
!     subroutine trsolv(a,b,c,f,x,lev0,lev1,k1,k2,lon0,lon1,lonmax,lat,
!    |  idebug)
!
!     call trsolv(p_coef,q_coef,r_coef,f_rhs,fcomp_out(:,lon0:lon1,lat),
!    |  lev0,lev1,lev0,lev1,lon0,lon1,nlonp4,lat,0)
!
      call trsolv(p_coef,q_coef,r_coef,f_rhs,fcomp_out(:,lon0:lon1,lat),
     |  lev0,lev1,lev0,lev1-1,lon0,lon1,nlonp4,lat,0)

!     if (idebug > 0) write(6,"('Minor second lat scan after ',
!    |  ' trsolv: lat=',i3)") lat
!
!     if (idebug > 0)
!    |  call addfsech('MNR_SOLV'  ,' ',' ',fcomp_out(:,lon0:lon1,lat),
!    |    lon0,lon1,nk,nk,lat)

      enddo ! lat=lat0,lat1
!----------------------- End second latitude scan ----------------------
!
! Wave filter the minor species output. Filter_minor does 3d gather/scatter
! for fft, so it is isolated from the latitude loop. Pass task subdomain,
! not including ghost cells.
!
!     if (idebug > 0) 
!    |  write(6,"('minor call filter_minor: nlat=',i4,' kut=',/,
!    |    (10i6))") nlat,kut

      call filter_minor(fcomp_out(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1,kut)
!
!----------------------- Begin third latitude scan --------------------
      do lat=lat0,lat1
!       write(6,"('Minor third latitude scan: lat=',i3)") lat

!       if (idebug > 0)
!    |    call addfsech('MNR_FILT'  ,' ',' ',fcomp_out(:,lon0:lon1,lat),
!    |    lon0,lon1,nk,nk,lat)
!
! Time smoothing:
        do i=lonbeg,lonend
          do k=lev0,lev1
            fcomp_tm1_out(k,i,lat) = dtsmooth*fcomp(k,i,lat)+
     |        dtsmooth_div2*(fcomp_tm1(k,i,lat)+fcomp_out(k,i,lat))
          enddo ! k=lev0,lev1
        enddo ! i=lonbeg,lonend

!       if (idebug > 0)
!    |    call addfsech('MNR_SMOO'  ,' ',' ',
!    |      fcomp_tm1_out(:,lon0:lon1,lat),lon0,lon1,nk,nk,lat)
!
      enddo ! lat=lat0,lat1
!----------------------- End third latitude scan ----------------------
!
! Periodic points for fcomp_out and fcomp_tm1_out:
#ifdef MPI
      call mp_periodic_f3d(fcomp_out(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1)
      call mp_periodic_f3d(fcomp_tm1_out(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1)
#endif
!
!----------------------- Begin fourth latitude scan --------------------
      do lat=lat0,lat1
!       if (idebug > 0) 
!    |    write(6,"('Minor fourth latitude scan: lat=',i3)") lat
!
! Insure density > 0.
! There is a much more elaborate check in earlier versions, which
! is not used here. This results in very small "diamond diffs" in
! N4S, and no diffs at all in NO (wrt tgcm15).
!
      do i=lon0,lon1
        do k=lev0,lev1
          if (fcomp_out(k,i,lat) < small) then
!           write(6,"('minor: fcomp_out < small: lat=',i2,' i=',
!    |        i2,' k=',i2,' fcomp_out=',e12.4)") 
!    |        lat,i,k,fcomp_out(k,i,lat)
            fcomp_out(k,i,lat) = small
          endif
          if (fcomp_tm1_out(k,i,lat) < small) then
!           write(6,"('minor: fcomp_tm1_out < small: lat=',i2,' i=',
!    |        i2,' k=',i2,' fcomp_out=',e12.4)") 
!    |        lat,i,k,fcomp_tm1_out(k,i,lat)
            fcomp_tm1_out(k,i,lat) = small
          endif
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1

!     if (idebug > 0) then
!       call addfsech('MNR_OUT'  ,' ',' ',
!    |    fcomp_out(:,lon0:lon1,lat),lon0,lon1,nk,nk,lat)
!       call addfsech('MNR_TM1'  ,' ',' ',
!    |    fcomp_tm1_out(:,lon0:lon1,lat),lon0,lon1,nk,nk,lat)
!     endif
!
! End fourth and final latitude loop:
      enddo ! lat=lat0,lat1
#ifdef VT
!     code = 114 ; state = 'minor' ; activity='ModelCode'
      call vtend(114,ier)
#endif
!
      end subroutine minor
!-----------------------------------------------------------------------
      subroutine filter_minor(fout,lev0,lev1,lon0,lon1,lat0,lat1,kut)
!
! Filter minor species (called from minor after trsolv)
!
      use params_module,only: nlonp4,nlon,nlat
      use filter_module,only: filter2
!
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
      integer,intent(in) :: kut(nlat)
      real,intent(inout) :: fout(lev0:lev1,lon0:lon1,lat0:lat1)
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
! Local:
      integer :: i,j,ier
      real :: f_ik(nlonp4,lev0:lev1),f_kij(lev0:lev1,nlonp4,lat0:lat1)
      real :: fmin,fmax
!
#ifdef VT
!     code = 123 ; state = 'filter_minor' ; activity='Filtering'
      call vtbegin(123,ier)
#endif
!     write(6,"('Enter filter_minor')")
!
! Define lons in f_ki from current task:
      f_kij = 0.
      do j=lat0,lat1
        do i=lon0,lon1
          f_kij(:,i,j) = fout(:,i,j)
        enddo
      enddo
!
#ifdef MPI
!
! Gather longitudes into tasks in first longitude column of task table
!   (leftmost of each j-row) for global fft. (i.e., tasks with mytidi==0
!   gather lons from other tasks in that row). This includes all latitudes.
!
      call mp_gatherlons_f3d(f_kij,lev0,lev1,lon0,lon1,lat0,lat1,1)
#endif
!
! Only leftmost tasks at each j-row of tasks does the global filtering:
      if (mytidi==0) then
!
! Loop through subdomain latitudes:
        latscan: do j=lat0,lat1
          if (kut(j) >= nlon/2) cycle latscan
!
! Define 2d array with all longitudes for filter2:
          do i=1,nlonp4
            f_ik(i,:) = f_kij(:,i,j)
          enddo ! i=1,nlonp4
!
! Do the filtering (requires fft):
          call filter2(f_ik,lev0,lev1,kut(j),j)
!
! Return filtered array to f_kij:
          do i=1,nlonp4
            f_kij(:,i,j) = f_ik(i,:)
          enddo ! i=1,nlonp4
        enddo latscan ! j=lat0,lat1
      endif ! mytidi==0
#ifdef MPI
!
! Now leftmost task at each j-row must redistribute filtered data
! back to other tasks in the j-row (mytidi>0,mytidj) (includes latitude):
!
      call mp_scatterlons_f3d(f_kij,lev0,lev1,lon0,lon1,lat0,lat1,1)
#endif
!
! Return filtered array to fout at current task longitudes:
      do j=lat0,lat1
        do i=lon0,lon1
          fout(:,i,j) = f_kij(:,i,j)
        enddo
      enddo ! j=lat0,lat1
!
#ifdef VT
!     code = 123 ; state = 'filter_minor' ; activity='Filtering'
      call vtend(123,ier)
#endif
      end subroutine filter_minor
