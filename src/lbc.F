!
      module lbc
!
! Set lower boundary conditions of t,u,v,z.
!
      use params_module,only: nlat,nlonp4,nlevp1,dz,zibot
      use cons_module,only: t0,pi,atm_amu,gask,grav,freq_semidi,
     |  re,dlamda,tgrad,cs,cor,tn,freq_ann,dt
      use init_module,only: istep
      use addfld_module,only: addfld
      implicit none
!    
! Total lower boundary conditions returned by this module.
! (dimensioned at full global grid, but defined at subdomains only)
! 
      real,dimension(nlonp4,nlat) :: t_lbc, u_lbc, v_lbc, z_lbc
!
      complex,parameter :: ci=(0.,1.), expta=1.
      complex ::
     |  zb(nlat),zb2(nlat),zba(nlat),
     |  tb(nlat),tb2(nlat),tba(nlat),
     |  ub(nlat),ub2(nlat),uba(nlat),
     |  vb(nlat),vb2(nlat),vba(nlat),
     |  bnd(nlonp4),bnd2(nlonp4),bnda(nlonp4)
!
      real :: b(nlonp4,2,2),fb(nlonp4,2)            ! for bndry_comp
!
! Planetary and Kelvin waves:
      real :: zbplanet(nlonp4,nlat)  ! planetary wave (sub bndry_planetary)
      real :: zbkelvin(nlonp4,nlat)  ! kelvin wave (sub bndry_kelvin)
!
      real,parameter ::      ! background tn at lower boundary
! 9/15/06 new vtgcm value from old vtgcm2 at 94 km
!    |  tbound = 171.0,      ! (94 km vtgcm)
! 8/02/07 new vtgcm value from Seiff et al., 1980
     |  tbound = 230.0       ! (72 km vtgcm)
!
      real ::
     |  ubound(nlat),   ! RSZ winds at lower boundary (lat dependent)
     |  vbound(nlat)    ! Meridional winds at lower boundary (lat dependent)
!
      contains
!-----------------------------------------------------------------------
      subroutine init_lbc
      use cons_module,only: cs,sn
!
! Set lower boundary Hough modes 
! (these may or may not be used, see sub tuvz_lbc)
! This is called once per run from tgcm.
!
! Local:
      integer :: j
      real ::
     |  sfac   = 0.,      ! sfac : Scale factor for RSZ flow in thermosphere
     |  uequ   = 9000.0,  ! U(RSZ) flow at equator at 94.0 km (cm/sec)
     |  usbr40 = 7000.0,  ! U(SBR40) SBR flow at 40LAT at 94.0 km (cm/sec)
     |  v30    = 1000.0   ! V(30) ZM 30-LAT meridional flow at 94.0 km (cm/sec)

      call bndry_diurnal
      call bndry_semidiurnal
      call bndry_comp
!
!  11/15/12 btf: Calculation of u,vbound moved here from sub init_cons (cons.F)
!    (at this time, u,vbound will be zero because sfac is zero)
!
!  Calculate ZM latitude dependent ubound(j) and vbound(j) for use at lower
!  boundary of vtgcm (cm/sec). 
!  -- RSZ flow prescription for sensitivity tests : 12/06/06
!     (SBR model or mid-latitude jet model)
!  -- Meridional flow prescription for sensitivity tests : 12/06/06
!     (50% of full meridional flow near 30 LAT in both hemispheres)
!  -- SFAC to describe the relative srength of RSZ flow in thermosphere
!     (sensitivity tests needed to examine impacts on winds, structure & AG)
!
      do j=1,nlat
!  ** TEST #1 -----(SBR, retrograde) ------------
        ubound(j) =  uequ*cs(j)*sfac
        vbound(j) =  2.0*v30*sn(j)*sfac
!  ** TEST #2 -----(MID-LAT JETS)----------------
!       ubound(j) =  (uequ*cs(j)+usbr40*cs(j)*abs(sn(j)))*sfac
!       vbound(j) =  2.0*v30*sn(j)*sfac
      enddo ! j=1,nlat
!     write(6,"('init_lbc: ubound=',/,(6e12.4))") ubound
!     write(6,"('init_lbc: vbound=',/,(6e12.4))") vbound

      end subroutine init_lbc
!-----------------------------------------------------------------------
      subroutine tuvz_lbc
!
! Update lower boundary subdomains of T,U,V,Z
! This is called every timestep from advance
!
      use input_module,only: oxvgcm_ncfile,step,tideann,planetary,kelvin
      use hist_module,only: modeltime
      use init_module,only: iter,iday
      use mpi_module,only: lon0,lon1,lat0,lat1
      use oxvgcm,only: t_oxvgcm,u_oxvgcm,v_oxvgcm,z_oxvgcm,oxvgcm_nlev,
     |  oxvgcm_bndry
!
! Local:
      integer :: i,j
      real :: rstep,dday
      complex :: expt_sdi, expt_di, expta
!
! Calculate exponentials
      rstep = float(step)
      expt_sdi   = cexp(ci*freq_semidi*rstep*iter)
      expt_di    = cexp(ci*.5*freq_semidi*rstep*iter)
      dday = float(iday)+amod(float(iter)*float(step),86400.)/86400.
!
! Set background lower boundaries (these may be added to or overriden later)
! (t,u,vbound are in cons.F)
!
      t_lbc(lon0:lon1,lat0:lat1) = tbound
      
      do j=lat0,lat1
        u_lbc(lon0:lon1,j) = ubound(j)
        v_lbc(lon0:lon1,j) = vbound(j)

        do i=lon0,lon1
          z_lbc(i,j) = real(zb(j)*bnd(i)*expt_sdi) ! see addiag.F
!
! Add in effect of (1,1) tidal component to lbc (not in vtgcm)
!         z_lbc(i,j) = z_lbc(i,j)+real(zb2(j)*bnd2(i)*expt_di)
        enddo
      enddo
!
! Override background with Oxford VGCM lower boundaries:  
! (oxvgcm_ncfile is namelist input path to oxvgcm data file)
! -- Offset of temperatures to better approximate real atmosphere
! -- (in readdata)
!
      if (len_trim(oxvgcm_ncfile) > 0) then
        call oxvgcm_bndry(modeltime)
        do j=lat0,lat1
          do i=lon0,lon1
            t_lbc(i,j) = t_oxvgcm(i,j,1)
            u_lbc(i,j) = u_oxvgcm(i,j,1)
            v_lbc(i,j) = v_oxvgcm(i,j,1)
            z_lbc(i,j) = z_oxvgcm(i,j,1)
          enddo
        enddo
        if (istep==1) write(6,"('tuvz_lbc: defined lbc of t,u,v,z ',
     |    'with Oxford VGCM (oxvgcm) data')")
      endif
!
! Add in effect of annual tidal component to lbc of Z:
! (tideann is a 0/1 namelist input, zero by default)
!
      if (tideann > 0) then
        call bndry_annual(dday)
        do j=lat0,lat1
          do i=lon0,lon1
            z_lbc(i,j) = z_lbc(i,j)+real(zba(j)*bnda(i)*expta)
          enddo
        enddo
        if (istep==1)write(6,"('tuvz_lbc: added annual tide to z_lbc')")
      endif
!
! Add planetary wave perturbations to lbc of Z:
! (planetary is a 0/1 namelist input, zero by default)
!
      if (planetary > 0) then
        call bndry_planetary  ! zbplanet
        do j=lat0,lat1
          do i=lon0,lon1
            z_lbc(i,j) = z_lbc(i,j)+zbplanet(i,j)
          enddo
        enddo
        if (istep==1) write(6,"('tuvz_lbc: added planetary to z_lbc')")
      endif
!
! Add Kelvin wave perturbations to lbc of Z:
! (kelvin is a 0/1 namelist input, zero by default)
!
      if (kelvin > 0) then
        call bndry_kelvin     ! zbkelvin
        do j=lat0,lat1
          do i=lon0,lon1
            z_lbc(i,j) = z_lbc(i,j)+zbkelvin(i,j)
          enddo
        enddo
        if (istep==1) write(6,"('tuvz_lbc: added kelvin to z_lbc')")
      endif

!     write(6,"('tuvz_lbc: t_lbc min,max=',2e12.4)") 
!    |  minval(t_lbc(lon0:lon1,lat0:lat1)),
!    |  maxval(t_lbc(lon0:lon1,lat0:lat1))
!     write(6,"('tuvz_lbc: u_lbc min,max=',2e12.4)") 
!    |  minval(u_lbc(lon0:lon1,lat0:lat1)),
!    |  maxval(u_lbc(lon0:lon1,lat0:lat1))
!     write(6,"('tuvz_lbc: v_lbc min,max=',2e12.4)") 
!    |  minval(v_lbc(lon0:lon1,lat0:lat1)),
!    |  maxval(v_lbc(lon0:lon1,lat0:lat1))
!     write(6,"('tuvz_lbc: z_lbc min,max=',2e12.4)") 
!    |  minval(z_lbc(lon0:lon1,lat0:lat1)),
!    |  maxval(z_lbc(lon0:lon1,lat0:lat1))
!
! Save to secondary histories:
        call addfld("T_LBC","t_lbc at end of sub tuvz_lbc","deg K",
     |    t_lbc(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
        call addfld("U_LBC","u_lbc at end of sub tuvz_lbc","cm/s",
     |    u_lbc(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
        call addfld("V_LBC","v_lbc at end of sub tuvz_lbc","cm/s",
     |    v_lbc(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
        call addfld("Z_LBC","z_lbc at end of sub tuvz_lbc","cm",
     |    z_lbc(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)

      end subroutine tuvz_lbc
!-----------------------------------------------------------------------
      subroutine bndry_semidiurnal
!
! Lower boundary conditions for semi-diurnal tide:
! 2/00: 1998 spherepack lib code (sphpac.f) replaces old lib 
!       alfpac.f for legendre polynomials and Hough functions.
! This routine calculates complex ZB,TB,UB,VB.
!
      use input_module,only: tide
!
! Local:
      integer,parameter :: nalf=19, malf=2
      real :: p(nlat,nalf,malf),hough(nlat,5,malf),cp(nalf/2+1)
      complex :: dzb(nlat)
      real :: B(5,19),RL(5),BHOUR(5),rlamda,xdot(19),ydot(19),
     |  ptscal,theta,ptjm(2*nlat+1)
      integer :: n,jm,l,lm1,m,mm1,j,ld,i,nm1
!
! External:
      real,external :: sddot ! in util.F
!
      COMPLEX ZEE(5),CL(5),EXPDLM
      DATA B/
     2 0.969152, 0.0     , 0.216046, 0.0     , 0.093838,
     3 0.0     , 0.909763, 0.0     , 0.342113, 0.0     ,
     4-0.245226, 0.0     , 0.798445, 0.0     , 0.421218,
     5 0.0     ,-0.408934, 0.0     , 0.645517, 0.0     ,
     6 0.024633, 0.0     ,-0.543993, 0.0     , 0.464159,
     7 0.0     , 0.071127, 0.0     ,-0.643189, 0.0     ,
     8-0.001292, 0.0     , 0.139613, 0.0     ,-0.699495,
     9 0.0     ,-0.006673, 0.0     , 0.225090, 0.0     ,
     1 0.000042, 0.0     ,-0.019654, 0.0     , 0.320141,
     1 0.0     , 0.000394, 0.0     ,-0.043345, 0.0     ,
     2-0.000001, 0.0     , 0.001772, 0.0     ,-0.079831,
     3 0.0     ,-0.000016, 0.0     , 0.005401, 0.0     ,
     4 0.0     , 0.0     ,-0.000112, 0.0     , 0.012932,
     5 0.0     , 0.0     , 0.0     ,-0.000476, 0.0     ,
     6 0.0     , 0.0     , 0.000005, 0.0     ,-0.001490,
     7 0.0     , 0.0     , 0.0     , 0.000031, 0.0     ,
     8 0.0     , 0.0     , 0.0     , 0.0     , 0.000129,
     9 0.0     , 0.0     , 0.0     ,-0.000002, 0.0     ,
     1 0.0     , 0.0     , 0.0     , 0.0     ,-0.000009/
!
!     DATA RL/7.8519E5, 3.6665E5, 2.1098E5, 1.3671E5, 0.9565E5/
!
! Equivalent depths for mars (2,2:2,3:2,4:2,5:2,6)
! All scaled from those of earth. actually use only first 3.
!
      DATA RL/5.4790E5, 2.5583E5, 1.4721E5, 0.9539E5, 0.66745E5/ ! mtgcm
!
      bhour = tide(6:10)
      if (all(tide==0.)) goto 13
!****** Commented out due to the moving LB tbound is dependent on Lat ******
!!    DO N=1,5
!!      ZEE(N)=TIDE(N)*CEXP(CI*pi*BHOUR(N)/6.)
!!      CL(N)=CSQRT(CMPLX(gask/(atm_amu*grav*RL(N))*
!!   |    (TBOUND*2./7.)-0.25))-0.5*CI ! mtgcm
!!    enddo
      JM=2*nlat+1
!
C     ****     SET UP HOUGH FUNCTIONS
! 
! Using new (1998) spherepack (sphpac.f):
      do n=2,nalf+1
        nm1 = n-1
        do m=2,malf+1
          mm1=m-1
          call alfk(n,m,cp)
          do j=1,jm
            theta = float(j-1)*pi/float(jm-1)
            call lfpt(n,m,theta,cp,ptscal)
            ptjm(j) = ptscal
          enddo
          do j=1,nlat
            p(j,nm1,mm1) = ptjm(2*(nlat+1-j))
          enddo
        enddo
        DO J=1,nlat
!         P(J,nm1,2)=SQRT(FLOAT(n*(n+1)-6))*P(J,nm1,2)-2.*TN(J)*
!    |               P(J,nm1,1)
          P(J,nm1,2)=-SQRT(FLOAT(n*(n+1)-6))*P(J,nm1,2)-2.*TN(J)*
     |               P(J,nm1,1)*(-1.) ! mtgcm
        enddo
      enddo
!
! Original code, using old alfpac (was alfpac.f):
!     DO L=2,20
!       LM1=L-1
!       DO M=2,3
!         MM1=M-1
!         CALL LFK(L,M,CP,W)
!         CALL LFP(L,M,JM,CP,PT)
!         DO J=1,nlat
!           P(J,LM1,MM1) = PT(2*(nlat+1-J))
!         enddo
!       enddo
!       DO J=1,nlat
!         P(J,LM1,2)=SQRT(FLOAT(L*(L+1)-6))*P(J,LM1,2)-2.*TN(J)*
!    |               P(J,LM1,1)
!       enddo
!     enddo
!
! util.F: real function sddot(n,x,y)
      DO L=1,5
        DO LD=1,2
          DO J=1,nlat
            xdot(:) = p(j,:,ld)
            ydot(:) = b(l,:)
            HOUGH(J,L,LD)=sddot(19,xdot,ydot)
          enddo
        enddo
      enddo
C     ****      GENERATE ZB, UB, VB AND TB
      DO 5 J=1,nlat
        TB(J)=0.
        ZB(J)=0.
        DZB(J)=0.
    5 CONTINUE
      DO L=1,5
        DO J=1,nlat
          ZB(J)=ZB(J)+ZEE(L)*HOUGH(J,L,1)
          DZB(J)=DZB(J)+ZEE(L)*HOUGH(J,L,2)
          TB(J)=TB(J)+CI*atm_amu*grav/gask*ZEE(L)*CL(L)*HOUGH(J,L,1)
        enddo
      enddo
      DO 7 J=1,nlat
        UB(J)=freq_semidi*re*(1.-(COR(J)/freq_semidi)**2)
        VB(J)=CI*grav*(DZB(J)-2.*COR(J)/(freq_semidi*CS(J))*ZB(J))/UB(J)
        UB(J)=grav*(COR(J)/freq_semidi*DZB(J)-2./CS(J)*ZB(J))/UB(J)
    7 CONTINUE
      GO TO 11
   13 CONTINUE
C     ****      ZERO BOUNDARY CONDITION
      DO 12 J=1,nlat
        ZB(J)=0.
        TB(J)=0.
        UB(J)=0.
        VB(J)=0.
   12 CONTINUE
   11 CONTINUE
C     ****      CALCULATE LONGITUDINAL STRUCTURE
!     RLAMDA = -2.*dlamda       ! tiegcm,timegcm
      RLAMDA = pi-2.*dlamda     ! mtgcm
      BND(1)=CEXP(CI*2.*RLAMDA)
      EXPDLM=CEXP(CI*2.*dlamda)
      DO 9 I=2,nlonp4
        BND(I)=BND(I-1)*EXPDLM
    9 CONTINUE
!
! bf 5/30/01: these compare favourably with mtgcm2.
!     write(6,"('bndry_semidiurnal: ')")
!     write(6,"('  TB = ',/,(6e12.4))") tb
!     write(6,"('  UB = ',/,(6e12.4))") ub
!     write(6,"('  VB = ',/,(6e12.4))") vb
!     write(6,"('  ZB = ',/,(6e12.4))") zb
!     write(6,"('  BND = ',/,(6e12.4))") bnd
!
      end subroutine bndry_semidiurnal
!-----------------------------------------------------------------------
      subroutine bndry_diurnal
!
C     ****     TIDAL BOUNDARY CONDITION FOR DIURNAL GRAVITATIONAL MODE
C     ****     (1,1)
! 2/00: 1998 spherepack lib code (sphpac.f) replaces old lib 
!       alfpac.f for legendre polynomials and Hough functions.
! This routine calculates complex ZB2,TB2,UB2,VB2.
!
      use input_module,only: tide2
!
! Local:
      integer,parameter :: nalf=19, malf=2
      real :: p(nlat,nalf,malf),hough(nlat,5,malf),cp(nalf/2+1)
      complex :: dzb(nlat)
      real :: B(1,19),RL (1),BHOUR2(1),rlamda,xdot(19),ydot(19),
     |  ptscal,theta,ptjm(2*nlat+1)
      integer :: l,m,j,n,jm,ld,i
      COMPLEX ZEE(1),CL(1),EXPDLM
!
      DATA B/
     2 0.282710,
     3 0.0     ,
     4-0.638229,
     5 0.0     ,
     6 0.620521,
     7 0.0     ,
     8-0.336408,
     9 0.0     ,
     1 0.117021,
     1 0.0     ,
     2-0.028332,
     3 0.0     ,
     4 0.005042,
     5 0.0     ,
     6-0.000686,
     7 0.0     ,
     8 0.000074,
     9 0.0     ,
     1-0.000006/
      DATA RL/0.6909E5/
      real,external :: sddot ! in util.F
!
      bhour2 = tide2(2)
      if (all(tide2==0.)) goto 13
      DO 1 N=1,1
        ZEE(N)=TIDE2(N)*CEXP(CI*pi*BHOUR2(N)/12.)
        CL(N)=CSQRT(CMPLX(gask/(atm_amu*grav*RL(N))*
     |  (T0(1)*2./7.+(T0(2)-T0(1))/dz)-.25))-.5*CI
    1 CONTINUE
      JM=2*nlat+1
!
C     ****     SET UP HOUGH FUNCTIONS
!
! Using new (1998) spherepack (sphpac.f):
      do n=1,19
        do m=1,2
          call alfk(n,m,cp)
          do j=1,jm
            theta = float(j-1)*pi/float(jm-1)
            call lfpt(n,m,theta,cp,ptscal)
            ptjm(j)=ptscal
          enddo
          do j=1,nlat
            p(j,n,m) = ptjm(2*(nlat+1-j))
          enddo
        enddo
        DO J=1,nlat
          P(J,n,2)=SQRT(FLOAT(n*(n+1)-2))*P(J,n,2)-TN(J)*P(J,n,1)
        enddo
      enddo
!
! Original code, using old alfpac (was alfpac.f):
!     DO L=1,19
!       DO M=1,2
!         CALL LFK(L,M,CP,W)
!         CALL LFP(L,M,JM,CP,PT)
!         DO J=1,nlat
!           P(J,L,M)=PT(2*(nlat+1-J))
!         enddo
!       enddo
!       DO J=1,nlat
!         P(J,L,2)=SQRT(FLOAT(L*(L+1)-2))*P(J,L,2)-TN(J)*P(J,L,1)
!       enddo
!     enddo
!
! util.F: real function sddot(n,x,y)
      DO L=1,1
        DO LD=1,2
          DO J=1,nlat
            xdot(:) = p(j,:,ld)
            ydot(:) = b(l,:)
            HOUGH(J,L,LD)=sddot(19,xdot,ydot)
          enddo
        enddo
      enddo
C     ****      GENERATE ZB2, UB2, VB2 AND TB2
      DO 5 J=1,nlat
      TB2(J)=0.
      ZB2(J)=0.
      DZB(J)=0.
    5 CONTINUE
      DO 6 L=1,1
      DO 6 J=1,nlat
      ZB2(J)=ZB2(J)+ZEE(L)*HOUGH(J,L,1)
      DZB(J)=DZB(J)+ZEE(L)*HOUGH(J,L,2)
      TB2(J)=TB2(J)+CI*atm_amu*grav/gask*ZEE(L)*CL(L)*HOUGH(J,L,1)
    6 CONTINUE
      DO 7 J=1,nlat
      UB2(J)=.5*freq_semidi*re*(1.-(COR(J)/(.5*freq_semidi))**2)
      VB2(J)=CI*grav*(DZB(J)-COR(J)/(.5*freq_semidi*CS(J))*ZB2(J))/
     |  UB2(J)
      UB2(J)=grav*(COR(J)/(.5*freq_semidi)*DZB(J)-1./CS(J)*ZB2(J))/
     |  UB2(J)
    7 CONTINUE
      GO TO 11
   13 CONTINUE
C     ****      ZERO BOUNDARY CONDITION
      DO 12 J=1,nlat
      ZB2(J)=0.
      TB2(J)=0.
      UB2(J)=0.
      VB2(J)=0.
   12 CONTINUE
   11 CONTINUE
C     ****      CALCULATE LONGITUDINAL STRUCTURE
      RLAMDA = -2.*dlamda
      BND2(1)=CEXP(CI*RLAMDA)
      EXPDLM=CEXP(CI*dlamda)
      DO 9 I=2,nlonp4
      BND2(I)=BND2(I-1)*EXPDLM
    9 CONTINUE
      end subroutine bndry_diurnal
!-----------------------------------------------------------------------
      subroutine bndry_annual(dday)
!
! Tidal boundary condition for annual mode
! For 1998 spherepack lib code (sphpac.f)
! (replaces old alfpac.f)
!
      use input_module,only: tideann
      use init_module,only: iday
!
! Args:
      real,intent(in) :: dday
!
! Local:
      real,parameter :: scale=1.
      integer,parameter :: nalf=24, malf=2
      real :: p(nlat,nalf,malf),hough(nlat,0:6,2),cp(nalf/2+1)
      complex :: dzb(nlat),zzb(nlat)
      real :: b(6,24),rl(0:6),rt2,rm,factor,xdot(24),ydot(24),
     |  pi,ptscal,ptjm(2*nlat+1),theta,w(nlat)
      integer :: n,l,m,i,k,jm,mp1,ld,j
      complex cc(0:6,0:6),cl(0:6),expt
!
      data ((b(i,j),i = 1,6),j = 1,12)/
     |-0.882922, 0.000000,-0.345087, 0.000000,-0.202228, 0.000000,
     | 0.000000,-0.930826, 0.000000,-0.301357, 0.000000,-0.152720,
     |-0.466226, 0.000000, 0.567457, 0.000000, 0.407114, 0.000000,
     | 0.000000,-0.362673, 0.000000, 0.694431, 0.000000, 0.438014,
     |-0.055436, 0.000000, 0.711847, 0.000000,-0.163050, 0.000000,
     | 0.000000,-0.044983, 0.000000, 0.625545, 0.000000,-0.325772,
     |-0.002909, 0.000000, 0.225723, 0.000000,-0.749160, 0.000000,
     | 0.000000,-0.002773, 0.000000, 0.186467, 0.000000,-0.723674,
     |-0.000086, 0.000000, 0.034940, 0.000000,-0.435919, 0.000000,
     | 0.000000,-0.000103, 0.000000, 0.029425, 0.000000,-0.379254,
     |-0.000002, 0.000000, 0.003267, 0.000000,-0.122687, 0.000000,
     | 0.000000,-0.000003, 0.000000, 0.002928, 0.000000,-0.104008/
      data ((b(i,j),i = 1,6),j = 13,24)/
     | 0.0     , 0.000000, 0.000206, 0.000000,-0.021267, 0.000000,
     | 0.0     , 0.0     , 0.000000, 0.000202, 0.000000,-0.018228,
     | 0.0     , 0.0     , 0.000009, 0.000000,-0.002540, 0.000000,
     | 0.0     , 0.0     , 0.000000, 0.000010, 0.000000,-0.002252,
     | 0.0     , 0.0     , 0.0     , 0.000000,-0.000223, 0.000000,
     | 0.0     , 0.0     , 0.0     , 0.0     , 0.000000,-0.000208,
     | 0.0     , 0.0     , 0.0     , 0.0     ,-0.000015, 0.000000,
     | 0.0     , 0.0     , 0.0     , 0.0     , 0.000000,-0.000015,
     | 0.0     , 0.0     , 0.0     , 0.0     ,-0.000001, 0.000000,
     | 0.0     , 0.0     , 0.0     , 0.0     , 0.000000,-0.000001,
     | 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.000000,
     | 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     /
      data (rl(n),n=1,6)/
     |  -10.8409E5,-7.0243E5,-2.4874E5,-1.9696E5,-1.0694E5,
     |  -0.9119E5/
!
! Annual boundary coefficients for lower boundary at z = -7.
!
      data((cc(k,n),n=0,6),k=0,3)/
     |  ( 0.136291E+03, 0.000000E+00),(-0.200536E-01, 0.000000E+00),
     |  ( 0.423456E+00, 0.000000E+00),(-0.143623E-02, 0.000000E+00),
     |  ( 0.262889E+00, 0.000000E+00),( 0.365122E-02, 0.000000E+00),
     |  ( 0.102716E+00, 0.000000E+00),( 0.826209E+00,-0.116536E-01),
     |  (-0.455993E+00,-0.311380E-02),( 0.271258E-01,-0.373420E-01),
     |  (-0.433095E-01, 0.165546E-01),( 0.713386E-02,-0.112538E-01),
     |  (-0.153849E-01, 0.103063E-01),( 0.214366E-02,-0.570878E-02),
     |  (-0.316298E+00, 0.115053E+00),(-0.159072E-01,-0.245495E-02),
     |  ( 0.302211E+00,-0.132446E-01),( 0.230750E-02, 0.170566E-03),
     |  ( 0.100434E+00,-0.299227E-02),( 0.264555E-02,-0.137723E-03),
     |  ( 0.499098E-01,-0.110255E-02),(-0.520584E-03,-0.114124E-02),
     |  ( 0.178599E-01, 0.561092E-02),( 0.557591E-03, 0.176165E-02),
     |  (-0.148151E-02, 0.749397E-03),(-0.617325E-03, 0.778743E-03),
     |  (-0.530835E-03, 0.641767E-03),(-0.964206E-03, 0.551394E-03)/
      data((cc(k,n),n=0,6),k=4,6)/
     |  (-0.137927E-02, 0.866386E-03),(-0.242825E-02, 0.441184E-03),
     |  ( 0.120715E-02,-0.136729E-02),( 0.122657E-03, 0.316213E-04),
     |  ( 0.390769E-03,-0.162978E-03),( 0.378377E-03,-0.195668E-04),
     |  ( 0.366912E-03,-0.681579E-04),(-0.470068E-03,-0.118650E-05),
     |  ( 0.120025E-02,-0.797459E-03),( 0.622700E-03,-0.424648E-04),
     |  (-0.537275E-03, 0.101658E-03),( 0.222407E-03,-0.828812E-05),
     |  (-0.209097E-03, 0.828365E-04),( 0.945940E-04, 0.317248E-04),
     |  ( 0.341903E-03, 0.192246E-04),( 0.129833E-03, 0.247156E-04),
     |  (-0.610206E-03, 0.591081E-06),(-0.102160E-03,-0.434110E-04),
     |  ( 0.196672E-04,-0.305687E-04),(-0.905354E-04,-0.813929E-04),
     |  ( 0.569460E-05,-0.116661E-03)/
      real,external :: sddot ! in util.F
!
      rt2 = sqrt(2.)
      if (tideann==0) then 
!    
! Zero boundary condition, except for geopotential:
!    
        do j = 1,nlat
          zba(j) = cc(0,0)*1./rt2*1.e5
          tba(j) = 0.
          uba(j) = 0.
          vba(j) = 0.
        enddo
        bnda = 1. ! whole array op
        return
      endif
      jm = 2*nlat+1
!    
! Height variation
!    
      do n = 1,6
        cl(n) = -csqrt(cmplx(gask/(atm_amu*grav*rl(n))*
     |    (t0(1)*2./7.+(t0(2)-t0(1))/dz)-.25))-.5*ci
      enddo
      cl(0) = 0.
!
! Set up legendre polynomials
!
! Using new (1998) spherepack (sphpac.f):
! (nalf=24, see cbndrya.h)
!
      pi = 4.*atan(1.)
      do n=1,nalf
        do m=0,1
          mp1=m+1
          call alfk(n,m,cp)
          do j=1,jm
            theta = float(j-1)*pi/float(jm-1)
            call lfpt(n,m,theta,cp,ptscal)
            ptjm(j) = ptscal
          enddo
          do j=1,nlat
            p(j,n,mp1) = ptjm(2*(nlat+1-j))
          enddo
        enddo
        m = 0
        rm = float(m)
        do j=1,nlat
          p(j,n,2)=sqrt(float(n*(n+1)-m*(m+1)))*p(j,n,2)-rm*tn(j)*
     |      p(j,n,1)
        enddo
      enddo
!
! Now evaluate hough functions
!
      do l=1,6
        do ld=1,2
          do j=1,nlat
            xdot(:) = p(j,:,ld)
            ydot(:) = b(l,:)
            hough(j,l,ld)=sddot(24,xdot,ydot)
          enddo
        enddo
      enddo
!
! Hough function of order zero
!
      do j=1,nlat
        hough(j,0,1) = 1./rt2
        hough(j,0,2) = 0.
      enddo
!
! Generate zba, tba, uba, vba
!
      do j=1,nlat
        zba(j) = 0.
        tba(j) = 0.
        uba(j) = 0.
        vba(j) = 0.
      enddo
!
! Summation over frequency, k, calculation of phase factor
!
      do k = 0,6
        expt = cexp(ci*float(k)*(dday-1.)*86400.*freq_ann)
!
! Summation over order, n
!
        do n = 0,6
          factor = scale
          if (k.eq.0.and.n.eq.0) factor = 1.
          do j = 1,nlat
            zzb(j) = cc(k,n)*hough(j,n,1)*expt*1.e5*factor
            dzb(j) = cc(k,n)*hough(j,n,2)*expt*1.e5*factor
            w(j) = grav/(re*((float(k)*freq_ann)**2-cor(j)**2))
            zba(j) = zba(j)+zzb(j)
            tba(j) = tba(j)+ci*atm_amu*grav/gask*cl(n)*zzb(j)
            uba(j) = uba(j)+w(j)*(cor(j)*dzb(j)-rm*float(k)*freq_ann/
     |        cs(j)*zzb(j))
            vba(j) = vba(j)+ci*w(j)*(float(k)*freq_ann*dzb(j)-rm*cor(j)/
     |        cs(j)*zzb(j))
          enddo
        enddo
      enddo
!
! Longitudinal structure
!
      bnda = 1. ! whole array op
      end subroutine bndry_annual
!-----------------------------------------------------------------------
      subroutine bndry_comp
C     ****
C     ****     CALCULATE MATRICES B(ZIMXP,2,2) AND VECTORS FB(ZIMXP,2)
C     ****       REPRESENTING THE LOWER BOUBNDARY CONDITION IN COMP,
C     ****       WHERE PSI1 AND PSI2 ARE CALCULATED:
C     ****
C     ****         PSI(K=-1/2) = B * PSI(K=1/2) + FB
C     ****
C     ****     BNDCMP CALLS THE SUBROUTINE BNDEF TO DEFINE THE 2 X 2
C     ****       MATRICES E, F AND THE 2 VECTOR G IN THE GENERAL
C     ****       LOWER BOUNDARY CONDITION:
C     ****
C     ****         E * D(PSI)/DS + F * PSI + G = 0.
C     ****
C     ****         WHERE:
C     ****           PSI = |PSI1| AND THE BOUNDARY CONDITION IS APPLIED
C     ****                 |    |
C     ****                 |PSI2|
C     ****
C     ****           AT LEVEL ZERO
C     ****
C     ****     THIS SUBROUTINE THEN EVALUATES B AND FB FROM:
C     ****
C     ****       B = (E/DS - F/2.)**(-1) * (E/DS + F/2.)
C     ****
C     ****       FB = (E/DS - F/2.)**(-1) * G
C     ****
!
! Local:
      real :: EE(nlonp4,2,2),FF(nlonp4,2,2),GG(nlonp4,2),WM1(nlonp4,2,2)
     |  ,WM2(nlonp4,2,2),WM3(nlonp4,2,2),WV1(nlonp4,2),WS1(nlonp4)
      integer :: l,m,i
C     ****
C     ****     CALL BNDEF TO DEFINE E, F AND G IN S1, S2 AND S3
C     ****
C     ****
      CALL BNDEF(ee,ff,gg)
C     ****
C     ****     WM1 = (E/DS - F/2.)
C     ****
C     ****     WM1 = (E/DS + F/2.)
C     ****
      DO 1 L = 1,2
        DO 1 M = 1,2
          DO 1 I = 1,nlonp4
            WM1(I,L,M) = EE(I,L,M)/dz-FF(I,L,M)/2.
            WM2(I,L,M) = EE(I,L,M)/dz+FF(I,L,M)/2.
    1 CONTINUE
C     ****
C     ****     WM3 = WM1**(-1)
C     ****
C     ****       WS1 = DET(WM1)
C     ****
      DO 2 I = 1,nlonp4
        WS1(I) = WM1(I,1,1)*WM1(I,2,2)-WM1(I,1,2)*WM1(I,2,1)
    2 CONTINUE
C     ****
C     ****     NOW INVERSE OF WM1 IN WM3
C     ****
      DO 3 I = 1,nlonp4
        WM3(I,1,1) = WM1(I,2,2)/WS1(I)
        WM3(I,1,2) = -WM1(I,1,2)/WS1(I)
        WM3(I,2,1) = -WM1(I,2,1)/WS1(I)
        WM3(I,2,2) = WM1(I,1,1)/WS1(I)
    3 CONTINUE
C     ****
C     ****     B = WM3 * WM2
C     ****
      DO 4 L = 1,2
        DO 4 M = 1,2
          DO 4 I = 1,nlonp4
            B(I,L,M) = WM3(I,L,1)*WM2(I,1,M)+WM3(I,L,2)*WM2(I,2,M)
    4 CONTINUE
C     ****
C     ****     FB = WM3 * G
C     ****
      DO 5 L = 1,2
        DO 5 I = 1,nlonp4
          FB(I,L) = WM3(I,L,1)*GG(I,1)+WM3(I,L,2)*GG(I,2)
    5 CONTINUE
      RETURN
      end subroutine bndry_comp
!-----------------------------------------------------------------------
      SUBROUTINE BNDEF(ee,ff,gg)
C     ****
C     ****     BNDEF DEFINES THE LOWER BOUNDARY CONDITION FOR THIS
C     ****       VERSION OF THE MODEL
C     ****
C     ****     THE LOWER BOUNDARY CONDITION FOR COMP IS:
C     ****
C     ****       E * D(PSI)/DS + F * PSI +G = 0.
C     ****
C     ****         WHERE:
C     ****           PSI = VECTOR(PSI1,PSI2)
C     ****           E AND F ARE 2 X 2 MATRICES
C     ****           G = VECTOR(G1,G2)
C     ****
C     ****           E, F AND G MAY BE FUNCTIONS OF LATITUDE & LONGITUDE
C     ****
C     ****       THIS SUBROUTINE DEFINES E, F AND G for BNDCMP
C     ****
!
! Args:
      real,intent(out) :: EE(nlonp4,2,2),FF(nlonp4,2,2),GG(nlonp4,2)
!
! Local:
!  *  Constant PSI LBC
!     -- Modify now with PCE solution for [O] and [CO] at 72 km.
!     -- PSIO:  based on Yung and DeMore (1999) solution at ~70 km.
!     -- PSICO: based on actual VTGCM PCE solution at zlb = -16.0
!     --      : close to Yung and DeMore (1999) PCE solution at ~58 km
!     real,parameter :: alfa1=2.0E-06, alfa2=1.0E-03   ! old
!     real,parameter :: alfa1=1.2E-10, alfa2=2.9E-05
!     real,parameter :: alfa1=1.37E-10, alfa2=2.64E-04
      real,parameter :: alfa1=1.37E-10, alfa2=2.9E-05
      integer :: i
C     ****
C     ****     IN TIGCM AND TIEGCM:
C     ****
C     ****       E = |0.  0.|
C     ****           |      |
C     ****           |0.  1.|
C     ****
C     ****       F = |1.  1.|
C     ****           |      |
C     ****           |0. -1.|
C     ****
C     ****       G = |-ALFA|
C     ****           |     |
C     ****           |   0.|
C     ****
C     ****       WHERE:
C     ****         ALFA = 0.22 + 0.014 = 0.234
C     ****
C     ****     IN VTGCM AND MMGCM: NO/LITTLE DATA AT L.BOUNDARY
C     ****
C     ****       E = |0.  0.|
C     ****           |      |
C     ****           |0.  0.|
C     ****
C     ****       F = |1.  0.|
C     ****           |      |
C     ****           |0.  1.|
C     ****
C     ****       G = |-ALFA1|
C     ****           |      |
C     ****           |-ALFA2|
C     ****
C     ****       WHERE GLOBAL MEAN DENSITY LBCS: MASS MIXING
C     ****         ALFA1 = 1.82E-04    PHIO
C     ****         ALFA2 = 1.27E-03    PHICO
C     ****
      DO 1 I = 1,nlonp4
	EE(I,1,1) = 0.
	EE(I,1,2) = 0.
	EE(I,2,1) = 0.
	EE(I,2,2) = 0.
	FF(I,1,1) = 1.
	FF(I,1,2) = 0.
	FF(I,2,1) = 0.
	FF(I,2,2) = 1.
C     **** GM SPECIFIED O & CO AFTER KRASNOPOLSKY [1993]
	GG(I,1) = -ALFA1
	GG(I,2) = -ALFA2
    1 CONTINUE
      RETURN
      end subroutine bndef
!-----------------------------------------------------------------------
      subroutine bndry_planetary
      use init_module,only: iter,glat
      use cons_module,only: pi,dlamda,dt
!
! Calculate contribution to ZB from planetary waves
!
! 07/31/09 btf: This version obtained from Hanli.
! 11/15/12: This routine was taken from timegcm
!
!     real,parameter :: amplan=0.0, t0plan=6.0E+5
      real,parameter :: amplan=5.00E4, t0plan=6.0E+5
!     real,parameter :: amplan=11.0E4, t0plan=2.5E+5
!     data amplan,t0plan/ 5.0E4,2.5E+5/
!     data amplan,t0plan/ 11.0E4,2.5E+5/
!
!     real,parameter :: frqplan=1.212e-5  ! 6 day
      real,parameter :: frqplan=3.356e-5   ! 52 hour
!
! Local:
      real :: time,fac1
      integer :: j,i,istartlat,iendlat,iterstart,iterend
      real,parameter :: startlat=32.5
      real,parameter :: endlat=82.5
      real,parameter :: lonphs=-140./180.*3.14159
!
!      iterstart = (150.*86400.+.1)/dt
!      iterstart = 0.0
!      iterstart = (15.*86400.+.1)/dt
!      iterend = (25.*86400.+.1)/dt     ! caseB
      iterstart = (30.*86400.+.1)/dt
      iterend = (40.*86400.+.1)/dt

      istartlat = 0
      iendlat = 0
      do j=1,nlat-1
        if (startlat >= glat(j) .and. startlat < glat(j+1)) istartlat=j
        if (endlat > glat(j) .and. endlat <= glat(j+1)) iendlat=j
      enddo
      if (istartlat==0 .or. iendlat==0) then
        write(6,"(/,'>>> bndry_planetary: could not find index to',
     |    ' startlat=',f8.2,': glat=',/,(6f8.2))") startlat,glat
        call shutdown('startlat')
      endif
      zbplanet = 0. ! init whole-array
!
! 10/19/10 btf and Jia Yue: Always define time (sub bndry_planetary)
      time = (iter-iterstart)*dt
      if (iter < iterstart .or. iter > iterend) then
         fac1 = exp(-(time/t0plan)**2)
      else
         fac1 = 1.
      endif
!
! 10/8/04 btf: To implement planetary waves, uncomment below code.
!
      do j = istartlat,iendlat
        do i = 1,nlonp4
! **** PLANETARY WAVE 1
          zbplanet(i,j) = amplan*(sin((glat(j)-30.)*pi/(60.-
     |           (glat(j)-60.)/3.)))**2*
     |           fac1*cos(-pi+(i-3)*dlamda-lonphs)   ! i=3 corresponds to -pi

!         zbplanet(i,j) = amplan*(sin((glat(j)-30.)*pi/(60.-
!    |           (glat(j)-60.)/3.)))**2*
!    |           (1.-exp(-time/t0plan))*sin(1.*(-pi+(i-3)*dlamda))

!         zbplanet(i,j) = amplan*(sin((glat(j)-30.)*pi/(60.-
!    |           (glat(j)-60.)/3.)))**2*
!    |           sin(1.*(-pi+(i-3)*dlamda))
!         zbplanet(i,j) = 0.
!
! **** PLANETARY WAVE 2
!         zbplanet(i,j) = amplan*(sin((glat(j)-30.)*pi/(60.-
!    |           (glat(j)-60.)/3.)))**2*
!    |           (1.-exp(-time/t0plan))*sin(2.*(-pi+(i-3)*dlamda))
!         zbplanet(i,j) = amplan*(sin((glat(j)-30.)*pi/(60.-
!    |           (glat(j)-60.)/3.)))**2*
!    |           sin(2.*(-pi+(i-3)*dlamda))
! **** PW 3
!          zbplanet(i,j) = amplan*(sin((glat(j)-30.)*pi/(60.-
!     |           (glat(j)-60.)/3.)))**2*
!     |           fac1*sin(-pi+3.*(i-3)*dlamda+frqplan*time)


        enddo
      enddo
!
! Top north latitude is zero:
      zbplanet(:,nlat) = 0.
      zbplanet(:,1) = 0.

!     write(6,"('bndry_planet: zbplanet min,max=',2e12.4)")
!    |  minval(zbplanet),maxval(zbplanet)

      end subroutine bndry_planetary
!-----------------------------------------------------------------------
      subroutine bndry_kelvin
!
! Calculate contribution to ZB from ultra-fast Kelvin Waves
!
      use init_module,only: iter,glat
      use cons_module,only: pi,dt
      use hist_module,only: nstep,modeltime ! for print only
!
! Local:
! 10/7/04 btf: amplan=1.0e5 resulted in NaN's in most fields after ~1/2 hour.
! 10/7/04 btf: amplan=1.0e4 resulted in NaN's in most fields after ~14 hours.
!  2/1/05 btf: changed t0plan from 2.5e4 to 2.5e5 for timegcm1.2
!
!     real,parameter :: amplan=1.0E5, t0plan=2.5E+4
!     real,parameter :: amplan=1.0E4, t0plan=2.5E+4 ! timegcm1.1
      real,parameter :: amplan=1.0E4, t0plan=2.5E+5 ! timegcm1.2
      real :: omg3d,wvx1,time
      integer :: j,i
!
      omg3d = 2.*pi/86400./3.    ! radian frequency of 3 days (rad/sec)
      wvx1 = 2.*pi/360.          ! wavenumber 1 (rad/degree)
      zbkelvin = 0.              ! init whole-array
      time = iter*dt

      do j = 2,nlat-1
        do i = 1,nlonp4
          zbkelvin(i,j) = zbkelvin(i,j)+
     |           amplan*exp(-(glat(j)/30.)**2)*
     |           sin(omg3d*time-wvx1*(i-3)*5.)*
     |           (1.-exp(-time/t0plan))
        enddo
!
!       write(6,"('bndry_kelvin: mtime=',3i4,' nstep=',i4,' dt=',f8.2,
!    |    ' iter=',i5,' time=',e12.4,' j=',i3,' zbkelvin(:,j)=',/,
!    |    (6e12.4))") modeltime(1:3),nstep,dt,iter,time,j,zbkelvin(:,j)
!
      enddo
      zbkelvin(:,1) = 0.    ! highest south latitude
      zbkelvin(:,nlat) = 0. ! highest north latitude

!     write(6,"('bndry_kelvin: zbkelvin min,max=',2e12.4)")
!    |  minval(zbkelvin),maxval(zbkelvin)

      end subroutine bndry_kelvin
!-----------------------------------------------------------------------
      end module lbc
