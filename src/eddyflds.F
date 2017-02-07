      module eddyflds_module
      implicit none
      contains
!-----------------------------------------------------------------------
      subroutine eddyflds(tn,barm,difk,dift,xmue,
     |  lev0,lev1,lon0,lon1,lat)
!
! ---------------------------------------------------------------------
! Eddy quantities for use association: comp.F, dt.F, duv.F, minor.F
! (1)	eddy diffusion coefficient (based upon von Zahn et al., 1980)
! (2)	eddy thermal conductivity (based on diffusion coeff and prndtl #) 
! (3)	eddy viscosity coefficient (based on diffusion coeff) 
! Modified last: 3/06/07 (SWB):  employ full von Zahn plus use cutoff
! Modified last: 10/25-29/07 (SWB):  specified SZA variation of difk
! (1) use idn(i) dy/nt index from chapman.F
! Modified last: 08/20/08 (SWB):  nightside Kzz after 1996 vtgcm2 scheme.
! Modified last: 10/31/13 (SWB):  use zhang et al (2012) scheme below 100 km
! Modified last: 04/16/15 (SWB):  day-night Kzz after Y&D(1982) 1-D code
!                              :  consistent with KINETICS/JPL, Parkinson (2015)
! ---------------------------------------------------------------------
!
      use params_module,only: nlevp1,nlonp4,zibot,dz
      use qrj_module,only: mk_xnmbari
      use chapman_module,only: idn_lonlat
      use fields_module,only: tlbc
      use cons_module,only: gask,grav,prndtl
      use addfld_module,only: addfld
      implicit none

! --------------------------------------------------------------------
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
!
! Input args: press vs longitude (2d (k,i) on midpoints):
!
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: 
     |  tn,    ! neutral temperature (deg K)
     |  barm   ! mean molecular weight
!
! Output args (eddy quantities (1/sec units):
! Outputs (for passing to comp, dt, duv, minor) : 
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) ::
     |  difk,    ! eddy diffusion coefficient (1/sec)
!    |  eddvk,   ! eddy diffusion coefficient (1/sec)
     |  dift,    ! eddy thermal conductivity coefficient (1/sec) 
     |  xmue     ! eddy viscosity coefficient (1/sec)

! --------------------------------------------------------------------
! Local:
      integer :: k,i,k0,k1,i0,i1,it
      real :: KMAXD, KMAXN, KMIN, SZP, ZPBK
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  xnmbari,    ! N*MBAR at interfaces
     |  xnmbar,    ! N*MBAR at midlevels 
     |  tni,        ! tn at interfaces
     |  xtot,       ! total number density at interfaces
!    |  eddvk,      ! eddy diffusion coefficient (1/sec)
     |  tauedd,     ! H**2/k - Eddy diffusion timescale (Earth Days)
!    |  taueddvz,    ! H**2/k - Eddy diffusion timescale (Earth Days)
     |  ht,         ! scale height at interfaces
     |  sht         ! pressure scale height at interfaces 
!
! For addfld calls:
      k0=lev0 ; k1=lev1 ; i0=lon0 ; i1=lon1
! For eddy diffusion coefficient parameters:
!     KMAXN= 2.0E+07; KMAXD=1.0E+06; KMIN = 1.0E+06; ZPBK = 0.0    ! Normal case
!     KMAXN= 6.5E+07; KMAXD=1.0E+06; KMIN = 3.0E+06; ZPBK = 0.0    ! Test 
!     KMAXN= 1.0E+07; KMAXD=1.0E+06; KMIN = 0.5E+06; ZPBK = 0.0    ! Minimal case
!     KMAXN= 1.0E+05; KMAXD=1.0E+06; KMIN = 1.0E+04; ZPBK = 0.0    ! Ext. Minimal case
!     KMAXN= 1.0E+05; KMAXD=1.0E+06; KMIN = 0.5E+06; ZPBK = 0.0    ! Ext. Minimal case #2
!     KMAXN= 4.0E+07; KMAXD=1.0E+06; KMIN = 2.0E+06; ZPBK = 0.0    ! Maximum case (best case)
!     KMAXN= 4.0E+07; KMAXD=1.0E+06; KMIN = 2.0E+06; ZPBK = 1.0    ! Maximum case
!     KMAXN= 4.0E+07; KMAXD=1.0E+06; KMIN = 2.0E+06; ZPBK = 2.0    ! Maximum case
!     KMAXN= 4.0E+07; KMAXD=1.0E+06; KMIN = 2.0E+06; ZPBK = -2.0   ! Maximum case
!     KMAXN= 8.0E+07; KMAXD=1.0E+06; KMIN = 2.3E+06; ZPBK = 0.0    ! Extreme max case
!     KMAXN= 1.0E+06; KMAXD=1.0E+06; KMIN = 1.0E+06; ZPBK = 0.0    ! Shutoff case
!     KMAXN= 1.0E+07; KMAXD=1.0E+06; KMIN = 2.3E+06; ZPBK = 0.0    ! 
!     KMAXN= 1.0E+07; KMAXD=1.0E+06; KMIN = 2.0E+06; ZPBK = 0.0    ! 
!     KMAXN= 1.0E+07; KMAXD=1.0E+06; KMIN = 3.5E+05; ZPBK = 0.0    ! 
!     KMAXN= 1.0E+07; KMAXD=1.0E+06; KMIN = 5.0E+04; ZPBK = 0.0    ! 
      KMAXN= 1.0E+07; KMAXD=1.0E+07; KMIN = 1.0E+05; ZPBK = 0.0    ! 1-D Y&D(1982)
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
! Calculate p0*e(-z)*barm/kT (1/(N*MBAR))at interfaces
!     subroutine mk_xnmbari(tni,barm,xnmbari,lev0,lev1,lon0,lon1,lat)
!
      call mk_xnmbari(tni,barm(:,lon0:lon1),
     |  xnmbari,
     |  lev0,lev1,lon0,lon1,lat)
! --------------------------------------------------------------------
!  Provide cm2/sec unit eddy diffusion coefficient (first pass)
!  Eddy coeff formulation from von Zahn et al., (1980) cutoff at levels
!  -- All on midlevels
      do i=lon0,lon1
        do k=lev0,lev1-1
          szp = zibot + .5*dz + (k-1)*dz
	  xtot(k,i) =  xnmbari(k,i)/barm(k,i)
! ---------------------------------------------------------------------
!    Day Variations Specified 
          if (idn_lonlat(i,lat)==1) then
!    von Zahn et al., (1980): MINIMAL DAYSIDE MIXING (<10% 1-D code)
!              difk(k,i)   = 0.7E+13/sqrt(xtot(k,i))
!              difk(k,i)   = 0.7E+12/sqrt(xtot(k,i))
!              if (difk(k,i) > 5.0E+06) difk(k,i) = 5.0E+06
!              if (difk(k,i) > 1.0E+06) difk(k,i) = 1.0E+06
!    Constant minimal dayside value for altitudes above 100 km
!              difk(k,i)   = KMAXD
!    von Zahn et al., (1980): MINIMAL DAYSIDE MIXING (For 1-D KINETICS)
               difk(k,i)   = 1.4E+13/sqrt(xtot(k,i))
               if (difk(k,i) > KMAXD) difk(k,i) = KMAXD
          else
! ---------------------------------------------------------------------
!    Night Variations Specified 
!    von Zahn et al., (1980): FULL NIGHTSIDE MIXING
!              eddvk(k,i)   = 1.4E+13/sqrt(xtot(k,i))
!              if (eddvk(k,i) > 8.0E+07) eddvk(k,i) = 8.0E+07
!    Old VTGCM2 scheme, with KMAXN at zpbk = 0.0 level
!              difk(k,i)  =  KMAXN*exp((SZP-ZPBK)/2.) 
!              if (difk(k,i) > KMAXN) difk(k,i) = KMAXN
!    von Zahn et al., (1980):  modified for  1-D KINETICS
               difk(k,i)   = 1.4E+13/sqrt(xtot(k,i))
               if (difk(k,i) > KMAXN) difk(k,i) = KMAXN
	  endif
! ---------------------------------------------------------------------
          if (difk(k,i) < KMIN) difk(k,i) = KMIN
!
!  Scale Height
          ht(k,i) = gask*tni(k,i)/(barm(k,i)*grav)
!
          tauedd (k,i) = (ht(k,i)**2.0/difk(k,i))/86400.

!         taueddvz (k,i) = (ht(k,i)**2.0/eddvk(k,i))/86400.
!     WRITE (*,*) 'taueddvz',taueddvz
!         tauedd (k,i) = (1./difk(k,i))/86400.
!
!
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
! 
!     call addfld('KZZ1',' ',' ',difk,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('KZZ1',' ',' ',difk(:,i0:i1),'lev',k0,k1,'lon',i0,
     | i1,lat)
      call addfld('TAUEDD',' ',' ',tauedd(:,i0:i1),'lev',k0,k1,
     | 'lon',i0,i1,lat)
!     call addfld('TAUEDV',' ',' ',taueddvz(:,lon0:lon1),
!    | 'lev',lev0,lev1,'lon',lon0,lon1,lat)
! --------------------------------------------------------------------
!  Provide 1/sec units for model usable quantities (second pass)
!  Divide eddy coefficient by pressure scale height squared
!  -- All on interfaces
      do i=lon0,lon1
        do k=lev0,lev1-1
	  sht(k,i) = gask*tni(k,i)/(barm(k,i)*grav)
          difk(k,i)   = difk(k,i)/sht(k,i)**2.0
          dift(k,i)   = difk(k,i)/prndtl
!         xmue(k,i)   = difk(k,i)
          xmue(k,i)   = difk(k,i)/prndtl
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
!     call addfld('KZZ2',' ',' ',difk,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('KZZ2',' ',' ',difk(:,i0:i1),'lev',k0,k1,'lon',i0,
     | i1,lat)
! --------------------------------------------------------------------
!
      end subroutine eddyflds
      end module eddyflds_module
! ---------------------------------------------------------------------
