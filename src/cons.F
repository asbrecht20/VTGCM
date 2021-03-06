!
#include <defs.h>
      module cons_module
      use params_module,only: dlat,dlon,dz,nlon,nlonp1,nlonp4,nlat,
     |  nlatp1,dlev,nlev,nmlat,nmlon,nmlonp1,zmbot,zmtop,zibot,zitop
      implicit none
!
! Define model constants. 
! -- Latest mods for vtgcm:  09/15/06  (Bougher)
! -- Latest mods for vtgcm:  11/03/06  (Bougher) prndtl = 10.
! -- Latest mods for vtgcm:  12/04/06  (Bougher) remove difk,dift,xmue
! -- Latest mods for vtgcm:  12/06/06  (Bougher) add ubound(j) and vbound(j)
! -- Latest mods for vtgcm:  01/17/07  (Bougher) sign for ubound(j) and vbound(j)
! -- Latest mods for vtgcm:  02/13/08  (Bougher) prndtl= 10, usr = 0.0, tbound = 230.0
!                                      gmlbc_oco = 0.; gmlbc_oco_dbl = 2*2.9E-05
!                                      grav = 850 at 70 km
! -- Latest mods for vtgcm:  05/15/08  (Bougher) psin2 = 0.069 mr
! -- Latest mods for vtgcm:  05/30/08  (Bougher) psin2 = 0.100 mr
! -- Latest mods for vtgcm:  02/18/16  (Bougher) psco = 3.67E-05 mmr
! -- Latest mods for vtgcm:  02/26/16  (Bougher) aeropk = max aerosol/UV-absorber layer heating
!                                      (units = K/day): 8K/day
! -- Latest mods for vtgcm:  03/03/16  (Bougher) psco = 1.88E-05 mmr
! -- Latest mods for vtgcm:  03/03/16  (Bougher) psco = 1.88E-05 mmr
! -- Latest mods for vtgcm:  03/23/16  (Bougher) aeropk = max aerosol/UV-absorber layer heating
!                                      (units = K/day): 40-50/day
! Parameter constants are cons_module module data and are accessed 
!   in subroutines via use-association. 
! Derived constants are cons_module module data, and are calculated
!   in sub init_cons (contained in cons_module). 
!   Sub init_cons is called by sub init (init_mod.F).
! Parameters referenced here are in params.h (s.a., dims.h)
!
! btf 6/13/13: Add mwt for so, so2
!
! Real parameter constants:
!
      real,parameter ::
     |  dzp  = dz,           ! alias for dz (also dlev)
     |  re   = 6.052e8,      ! venus radius (cm)                  C(51)
     |  re_inv = 1./re,      ! inverse of earth radius            C(52)
     |  avo  = 6.023e23,     ! avogadro number                    C(85)
     |  boltz = 1.38E-16,    ! boltzman's constant                C(84)
!    |  p0   = 5.0e-4,       ! standard pressure                  C(81)
     |  p0   = 5.0e-3,       ! standrard pressure for vtgcm       C(81)
     |  gask = 8.314e7,      ! gas constant                       C(57)
     |  secperday = 2.0997e+7,   ! secs per venus day
     |  secperhr = secperday/24.,! secs per venus hour (8.75e+5)
!
! dipmin should be same as sin10 (see magfield.F):
!
#if (NLAT==36 && NLON==72)
     |  dipmin = 0.17,       ! minimum mag dip angle (5.0 deg horizontal res)
#elif (NLAT==72 && NLON==144)
     |  dipmin = 0.24,       ! minimum mag dip angle (2.5 deg horizontal res)
#elif (NLAT==18 && NLON==72)              /* vtgcm single hemisphere */
     |  dipmin = 0.17,       ! minimum mag dip angle (2.5 deg horizontal res)
#else
      UNKNOWN NLAT,NLON      ! compilation will stop here if unknown res
#endif
     |  brn2d = 0.6,         !
     |  atm_amu = 44.0,      ! mean mass of venus surface atmosphere    C(24)
     |  shapiro = 3.0E-2,    ! shapiro smoother constant          C(26)
     |  dtsmooth = 0.95,     ! time smoothing constant            C(30)
     |  dtsmooth_div2 = 0.5*(1.-dtsmooth), !                      C(31)
     |  tgrad = 6.,          ! TN gradient from old /RUNMDS/ (tiegcm)
     |  nob(nlat) = 4.e6,    ! N(NO) LBC from old /RUNMDS/ 
     |  avto = 4.0e-12,      ! from old /RUNMDS/ (not used)
     |  hor(nlat) = .25,     ! horizontal variation of eddy diffusion and
                             !  eddy thrmal conductivity. If unity, value
                             !  of KE at bottom is 5.e-6 (old /RUNMDS/)
!    |  prndtl = 1.,         ! prandtl number
     |  prndtl = 10.,        ! prandtl number (reduce role of econd for Venus)
     |  evergs = 1.602e-12,  ! 1 eV = 1.602e-12 ergs
!    |  aeropk = 70.0,       ! (K/day) peak aerosol heating from Titov (95) and Parkinson (15)
!    |  aeropk = 60.0,       ! (K/day) peak aerosol heating from Titov (95) and Parkinson (15)
     |  aeropk = 50.0,       ! (K/day) peak aerosol heating from Titov (95) and Parkinson (15)
     |  tsurplus=5.11*evergs, ! surplus heat per event (ergs)    C(45)
     |  amu = 1.66e-24,      ! atomic mass unit (g)
     |  gmlbc_oco    = 0.0,  ! vtgcm-ext (O-mass mixing ratio at z-16) C(46)
     |  gmlbc_oco_dbl= 1.88E-05, ! vtgcm-ext (CO mass mixing ratio z-16) C(47): 1-D Kinetics
!    |  gmlbc_oco_dbl= 2.*1.88E-05, ! vtgcm-ext (CO mass mixing ratio z-16) C(47):  1-D Kinetics
!    |  psin2 = 0.01684,     ! global mean lbc of N2 (mixing ratio)
!    |  psin2 = 0.043,       ! global mean lbc of N2 (mixing ratio)
!    |  psin2 = 0.069,       ! global mean lbc of N2 (mixing ratio)
     |  psin2 = 0.100        ! global mean lbc of N2 (mixing ratio)
!
      integer :: nlonper=nlonp4 ! nlon + periodic points (alias for nlonp4)
!
! Many expressions require x/rmass, but its more efficient on some
! platforms to multiply rather than divide, so set rmassinv = 1./rmass 
! here, and use x*rmassinv in the code.
!
      real,parameter :: 
     |  rmass(3) = (/16.,28.,44./), ! o1,co,co2
     |  rmass_o2 = 32., rmass_o1  = 16., rmass_n2  = 28., 
     |  rmass_o3 = 48., rmass_n4s = 14., rmass_n2d = 14., 
     |  rmass_no = 30., rmass_op  = 16., rmass_co2 = 44.,
     |  rmass_co = 28., rmass_so  = 48., rmass_so2 = 64.,
     |  rmass_so3 = 80.
      real,parameter ::
     |  rmassinv_o2  = 1./rmass_o2,
     |  rmassinv_o1  = 1./rmass_o1,
     |  rmassinv_n2  = 1./rmass_n2,
     |  rmassinv_o3  = 1./rmass_o3,
     |  rmassinv_n4s = 1./rmass_n4s,
     |  rmassinv_n2d = 1./rmass_n2d,
     |  rmassinv_no  = 1./rmass_no,
     |  rmassinv_op  = 1./rmass_op,
     |  rmassinv_co2 = 1./rmass_co2,
     |  rmassinv_co  = 1./rmass_co,
     |  rmassinv_so  = 1./rmass_so,
     |  rmassinv_so2 = 1./rmass_so2,
     |  rmassinv_so3 = 1./rmass_so3
!
! 2/00: these were in modsrc.snoe (tgcm13mt), but were unused.
! Low-energy protons:
!     real,parameter ::
!       alfalp = 10.,
!       efluxlp = 1.e-20
!
! Model derived constants (see sub init_cons in this module):
!
      real ::
     |  pi,             ! set with 4*atan(1)    C(110)
     |  rtd,            ! radians-to-degrees (180./pi)
     |  dtr,            ! degrees-to-radians (pi/180.)
     |  dphi,           ! delta lat (pi/nlat)   C(2)
     |  dphi_2div3,     ! 2./(3.*dphi)          C(12)
     |  dphi_1div12,    ! 1./(12.*dphi)         C(13)
     |  dlamda,         ! delta lon (2pi/nlon)  C(1)
     |  dlamda_2div3,   ! 2./(3.*dlamda)        C(10)
     |  dlamda_1div12,  ! 1./(12.*dlamda)       C(11)
     |  dt,             ! time step (secs)      C(4)
     |  dtx2,           ! 2*dt                  C(6)
     |  dtx2inv,        ! 1./(2*dt)             C(7)
     |  freq_3m3,       ! frequency of 2-day wave (rad/sec)       C(21)
     |  freq_semidi,    ! frequency of semidiurnal tide (rad/sec) C(23)
     |  freq_ann,       ! frequency of annual tide                C(25)
     |  expz(nlev+1),   ! exp(-z) at midpoints
     |  expzmid,        ! exp(-.5*dz)                             C(86)
     |  expzmid_inv,    ! 1./expzmid                              C(87)
     |  rmassinv(3),    ! inverse of rmass
     |  t0(nlev+1),     ! set by sub lowbound (bndry_mod.F)
     |  racs(nlat),     ! 1./(re*cs(lat))
     |  cs(-1:nlat+2),  ! cos(phi)
     |  sn(nlat),       ! sin(phi)
     |  tn(nlat),       ! tan(phi)
!    |  ubound(nlat),   ! RSZ winds at lower boundary (lat dependent)
!    |  vbound(nlat),   ! Meridional winds at lower boundary (lat dependent)
     |  cor(nlat),
     |  grav,           ! accel due to gravity (dependent on lower boundary)
     |  dzgrav          ! grav/gask C(65)
!
! 6/8/06 btf: difk and dift are calculated in mgw.F, and should be use-assoc
!             from there (see comp.F, minor.F, dt.F). xmue should be renamed
!             difv. After Ray calcs difv, all 3 can be removed from here.
! 12/04/06 swb: difk, dift, and xmue are calculated in eddyflds.F, & are use-assoc
!             from there (see comp.F, minor.F, dt.F, duv.F). 
!    |  difk(nlev+1),   ! background eddy diffusion
!    |  dift(nlev+1),   ! background thermal conductivity
!    |  xmue(nlev+1)    ! eddy viscosity (?)
!
! Constants for dynamo and electric field calculations:
! 5/10/06 btf: shift reference ht up from 80 to 90 km (s.a., kbotdyn)
      real,parameter :: h0 =8.0e6, r0 =re+h0 ! use mean earth radius
!     real,parameter :: h0 =9.0e6, r0 =re+h0 ! use 90km when kbotdyn ~= zp-8
      real :: 
     |  dlatg, dlong, dlatm, dlonm,
     |  ylatm(nmlat),    ! magnetic latitudes (radians)
     |  ylonm(nmlonp1),  ! magnetic longitudes (radians)
     |  rcos0s(nmlat),   ! cos(theta0)/cos(thetas)
     |  dt0dts(nmlat),   ! d(theta0)/d(thetas)
     |  dt1dts(nmlat),   ! dt0dts/abs(sinim) (non-zero at equator)
     |  table(91,2)
!
! Geographic grid in radians:
      real ::
     |  ylatg(0:nlatp1), ! geographic latitudes (radians)
     |  ylong(nlonp1)    ! geographic longitudes (radians)
!
! Critical colatitude limits for use of Heelis potential in dynamo:
      real,parameter :: 
     |  crit(2) = (/0.261799387, 0.523598775/)
!
! Kut is used in filtering longitudinal waves (see filter.F):  
!
#if (NLAT==36 && NLON==72)                 /* 5.0 deg horizontal resolution */
      integer,parameter :: kut(nlat) =
     |  (/1,2,3,5,6,7,9,10,11,13,14,15,17,17,17,17,17,17,17,17,17,17,17,
     |    17,15,14,13,11,10,9,7,6,5,3,2,1/)
!     integer :: kut(nlat) ! see sub set_wave_filter
#elif (NLAT==72 && NLON==144)              /* 2.5 deg horizontal resolution */
       integer,parameter :: kut(nlat) =
     |   (/1  ,1  ,2  ,2  ,4  ,4  ,8  ,8  ,10 ,10 ,12 ,12,
     |     15 ,15 ,18 ,18 ,22 ,22 ,26 ,26 ,30 ,30 ,32 ,32,
     |     34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34,
     |     34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34,
     |     32 ,32 ,30 ,30 ,26 ,26 ,22 ,22 ,18 ,18 ,15 ,15,
     |     12 ,12 ,10 ,10 ,8  ,8  ,4  ,4  ,2  ,2  ,1  ,1/)
#elif (NLAT==18 && NLON==72)              /* vtgcm single hemisphere */
       integer,parameter :: kut(nlat) =
     |   (/36,36,36,36,36,17,15,14,13,11,10,9,7,6,5,3,2,1/)
#else
      UNKNOWN NLAT,NLON      ! compilation will stop here if unknown res
#endif
!
       integer,parameter :: kut_5(36) =
     |  (/1,2,4,8,10,12,15,18,22,26,30,32,34,34,34,34,34,34,34,34,34,34,
     |    34,34,32,30,26,22,18,15,12,10,8,4,2,1/)
 
! If check_exp is set true, certain routines will use expo() (util.F)
! instead of exp(). expo checks for out of range arguments to the
! exponential, substituting large or small results if the argument
! is out of range. This avoids NaNS fpe's, but degrades performance.
! It will also produce slightly different results.
!
#ifdef DEBUG
      logical,parameter :: check_exp = .true.
#else
      logical,parameter :: check_exp = .false.
#endif
!
! Special pi for mag field calculations. If pi=4.*atan(1.) and code is
! linked with -lmass lib, then the last 2 digits (16th and 17th) of pi
! are different (56 instead of 12), resulting in theta0(j=49)==0., which 
! is wrong (should be .1110e-15).
!
      real,parameter ::
     |  pi_dyn=3.14159265358979312
!
      contains
!-----------------------------------------------------------------------
      subroutine init_cons
      use input_module,only: step
      use params_module,only: tgcm_version
!
! Set derived constants (this is called from sub init in init_module)
! -- Latest vtgcm mods:  09/15/06 (Bougher)
!
! Local:
      real :: expdz,phi
      real :: omega = -3.00E-7
!     real :: omega = 0.00
      integer :: k,i,j,js
!
      pi = 4.*atan(1.)                ! C(110)

      write(6,"('init_cons: pi=',e12.4)") pi

      rtd = 180./pi                   ! radians to degrees
      dtr = pi/180.                   ! degrees to radians
      dphi = pi/float(nlat)           ! C(2)
      dphi_2div3 = 2./(3.*dphi)       ! C(12)
      dphi_1div12 = 1./(12.*dphi)     ! C(13)
      dlamda = 2.*pi/float(nlon)      ! C(1)
      dlamda_2div3 = 2./(3.*dlamda)   ! C(10)
      dlamda_1div12 = 1./(12.*dlamda) ! C(11)
      t0(:) = 0.
!
! ------------------------------------------------------------------------ 
! Old way of computing difk, dift, xmue
! Eddy parameterizations now done in eddyflds.F (12/04/06) 
!
! ------------------------------------------------------------------------ 
!  --- Earth
!     difk(1) = 5.0e-6
!     dift(1) = 5.0e-6/prndtl
!     difk(1) = 5.0e-6
!  --- Venus (Kzz = 1.0E+07 cm2/sec, H = 3.8 km)
!     difk(1) = 6.9e-5
!     dift(1) = 6.9e-5/prndtl
!     xmue(1) = 6.9e-5
!     do k=2,nlev
!       difk(k) = difk(k-1)*expdz
!       xmue(k) = difk(k)
!       dift(k) = dift(k-1)*expdz
!     enddo
!     difk(nlev+1) = difk(nlev)*expdz
!     dift(nlev+1) = dift(nlev)*expdz
!     xmue(nlev+1) = difk(nlev)
! ------------------------------------------------------------------------ 
!
! expz(nlevp1) is exp(-zp) at midpoints (expz(nlev+1) not used):
!
      expz(:) = 0. ! init e(-z) at midpoints
      expz(1) = exp(-zmbot) ! zmbot is bottom midpoint level (see params.F)
      expdz = exp(-dlev)
      do k=2,nlev
        expz(k) = expz(k-1)*expdz
      enddo
      expzmid = exp(-.5*dlev)
      expzmid_inv = 1./expzmid
! ------------------------------------------------------------------------ 
      do i=1,3
        rmassinv(i) = 1./rmass(i)
      enddo
      js=-(nlat/2)
      do j=1,nlat
        phi=(j+js-.5)*dphi
        cs(j)=cos(phi)
        sn(j)=sin(phi)
        tn(j)=tan(phi)
        cor(j)=2.*omega*sn(j)
        racs(j) = 1./(re*cs(j))
      enddo
! ------------------------------------------------------------------------ 
!
! 11/15/12 btf: Lower boundary constants t,u,vbound are now in subroutine 
!               init_lbc in lbc.F
!
! ------------------------------------------------------------------------ 
! cs at 0, -1, nlat+1, and nlat+2 replace the old cssp and csnp:
      cs(-1) = -cs(2)
      cs(0)  = -cs(1)
      cs(nlat+1) = -cs(nlat)
      cs(nlat+2) = -cs(nlat-1)
      dt = float(step)                   ! was C(4)
      dtx2 = 2.*dt                       ! was C(6)
      dtx2inv = 1./dtx2                  ! was C(7)
      freq_3m3 = 2.*pi/(49.7789*60.*60.) ! was C(21)
      freq_semidi = 4.*pi/(24.*60.*60.)  ! was C(23)
      freq_ann = freq_semidi/(2.*365.25) ! was C(25)
!
! Set gravity according to lower boundary:
! Perform all plots using height = 0 flag (tgcmproc_f90)
!     grav = 835. ! (at 125 km in vtgcm; average of height range)
!     grav = 840. ! (at 94 km in vtgcm; z= -10 lower boundary)
      grav = 850. ! (at 70 km in vtgcm; z = -16 lower boundary)
      dzgrav = grav/gask   ! C(65)
!
!
! Set dynamo constants:
      call consdyn
!
! Report to stdout:
      write(6,"(  'Model = ',a)") tgcm_version
      write(6,"(/,'Set constants:')")
      write(6,"('  nlat=',i3,' nlon=',i3,' nlev=',i3)") nlat,nlon,nlev
      write(6,"('  dz = ',f5.2)") dz
      write(6,"('  zmbot, zmtop = ',2f8.3,
     |  ' (bottom,top midpoint levels)')") zmbot,zmtop
      write(6,"('  zibot, zitop = ',2f8.3,
     |  ' (bottom,top interface levels)')") zibot,zitop
      write(6,"('  dt = ',f8.2,' secs')") dt
      write(6,"('  grav = ',f10.2)") grav
      write(6,"('  freq_3m3 = ',e12.4,'  freq_semidi=',e12.4,
     |  '  freq_ann=',e12.4)") freq_3m3,freq_semidi,freq_ann
      write(6,"('  dipmin = ',f8.3)") dipmin
      write(6,"('  check_exp = ',l1)") check_exp
      write(6,"('  dlat=',f6.2,' dlon=',f6.2)") dlat,dlon
      write(6,"('  kut(nlat)=',/,(12i4))") kut

      end subroutine init_cons
!-----------------------------------------------------------------------
      subroutine consdyn
      use input_module,only: dynamo
!
! Set derived constants used in dynamo.
!
! Local:
      integer :: j,i,n
      real,parameter :: e=1.e-6, r1=1.06e7, alfa=1.668
      real :: 
     |  tanth0(nmlat),
     |  tanths(nmlat),
     |  theta0(nmlat),
     |  hamh0(nmlat)
      real :: dtheta,table2(91,3:5),tanths2
      real :: fac,rmin,rmax,rmag
!
! Set grid deltas:
      dlatg = pi/float(nlat)
      dlong = 2.*pi/float(nlon)
      dlatm = pi_dyn/float(nmlat-1) ! note use of pi_dyn
      dlonm = 2.*pi_dyn/float(nmlon) 
!
! Set geographic latitude array ylatg:
      do j=1,nlat
        ylatg(j) = -0.5*(pi-dlatg)+float(j-1)*dlatg
      enddo ! j=1,nlat
      ylatg(0) = -pi/2.+e
      ylatg(nlatp1) = pi/2.-e
!
! Set geographic longitude array ylong:
      do i=1,nlonp1
        ylong(i) = -pi+float(i-1)*dlong
      enddo ! i=1,nmlonp1
!
! Set magnetic latitudes ylatm and magnetic longitudes ylonm:
!
! ylatm is equally spaced in theta0, but holds corresponding value
!   of thetas.
      do j=1,nmlat
        theta0(j) = -pi_dyn/2.+float(j-1)*dlatm ! note use of pi_dyn
      enddo ! j=1,nmlat
      do j=2,nmlat-1
        tanth0(j) = abs(tan(theta0(j)))
        hamh0(j) = r1*tanth0(j)+r0*tanth0(j)**(2.+2.*alfa)/
     |    (1.+tanth0(j)**2)**alfa
        tanths(j) = sqrt(hamh0(j)/r0)
        ylatm(j) = sign(atan(tanths(j)),theta0(j))
        rcos0s(j) = sqrt((1.+tanths(j)**2)/(1.+tanth0(j)**2))
!
! If dynamo <= 0 -> no dynamo
! If dynamo >= 1 -> dynamo is called
!
        if (dynamo > 0) then
          tanths2  = tanths(j)**2
          dt1dts(j) = 
     |      (r0*sqrt(1.+4.*tanths2)*(1.+tanths2))/
     |      (r1*(1.+tanth0(j)**2)+2.*r0*tanth0(j)**(2.*alfa+1.)*
     |      (1.+alfa+tanth0(j)**2)/(1.+tanth0(j)**2)**alfa)
          dt0dts(j) = dt1dts(j)*2.*tanths(j)/sqrt(1.+4.*tanths2)
        endif
      enddo ! j=2,nmlat-1
!
! Magnetic poles:
      ylatm(1) = theta0(1)
      ylatm(nmlat) = theta0(nmlat)
      rcos0s(1) = 1.
      rcos0s(nmlat) = 1.
      dt0dts(1) = 1.
      dt0dts(nmlat) = 1.
!
! Magnetic longitudes:
      do i=1,nmlonp1
        ylonm(i) = -pi+float(i-1)*dlonm
      enddo ! i=1,nmlonp1
      dtheta = pi/(2.*90.)
!
      table(1,1) = 0. 
      table(1,2) = 0. 
      do i=2,91
        table(i,1) = table(i-1,1)+dtheta
      enddo
      do i=2,90
        table2(i,4) = tan(table(i,1))
        table(i,2) = table(i,1)
      enddo ! i=2,90
      table(91,2) = table(91,1)
!     table(91,2) = pi/2.
      do n=1,7
        do i=2,90
          table2(i,3) = table(i,2)
          table(i,2) = tan(table2(i,3))
          table2(i,5) = sqrt(r1/r0*table(i,2)+table(i,2)**(2.*(1.+alfa))
     |      /(1.+table(i,2)**2)**alfa)
          table(i,2) = table2(i,3)-(table2(i,5)-table2(i,4))*2.*
     |      table2(i,5)/(r1/r0*(1.+table(i,2)**2)+2.*table(i,2)**
     |      (2.*alfa+1.)*(1.+alfa+table(i,2)**2)/(1.+table(i,2)**2)**
     |      alfa)
        enddo ! i=2,90
      enddo ! n=1,7

!     write(6,"(/,'consdyn: ylatg =',/,(6e12.4))") ylatg
!     write(6,"(  'consdyn: ylong =',/,(6e12.4))") ylong
!     write(6,"(  'consdyn: ylatm =',/,(6e12.4))") ylatm
!     write(6,"(  'consdyn: rcos0s=',/,(6e12.4))") rcos0s
!     write(6,"(  'consdyn: dt0dts=',/,(6e12.4))") dt0dts
!
      end subroutine consdyn
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      end module cons_module
