!
! defs.h contains c-preprocessor macros.
!
#include <defs.h>
!
      module params_module
!
! Geographic grid parameters:
      integer,parameter :: 
     |  nlat = NLAT,                 ! number of latitudes
     |  nlon = NLON,                 ! number of longitudes
     |  nlev = NLEV,                 ! number of midpoint levels 
     |  nilev = NLEV,                ! number of interface levels
     |  nlonp4 = nlon+4,             ! includes periodic points in lon
     |  nlevp1 = nlev+1,
     |  nilevp1 = nilev+1,
     |  nlonp1=nlon+1, nlatp1=nlat+1,
     |  nlonp2=nlon+2, nlatp2=nlat+2
      real,parameter ::
     |  glat1 = GLAT1,               ! first latitude
     |  dlat  = DLAT,                ! delta latitude
     |  glon1 = GLON1,               ! first longitude
     |  dlon  = DLON,                ! delta longitude
     |  dlev  = (ZITOP-ZIBOT)/NLEV,  ! delta level (same for midpoints and interfaces)
     |  dz    = dlev,                ! delta level (alias for dlev)
     |  zmtop = ZMTOP,               ! top midpoint level
     |  zmbot = ZMBOT,               ! bottom midpoint level
     |  zitop = ZITOP,               ! top interface level
     |  zibot = ZIBOT                ! bottom interface level
      real,parameter ::
     |  spval = 1.e36
      integer,parameter ::
     |  ispval = 999
!
!     real,parameter :: 
! Venus: Rotation is about 12 deg per 5 earth days, or 0.1 deg per hour
!    |  secperday=1.296e+7,         ! seconds/day
!    |  secperhr =secperday/24.,    ! seconds/hr (540,000)
!    |  degperhr = 0.1              ! degrees rotation per hour
!
! Earth:
!    |  secperday=24.*60.*60.,      ! seconds/day
!    |  secperhr =60.*60.,          ! seconds/hr
!    |  degperhr =15.               ! rotation/hr
!
! Magnetic grid:
      integer,parameter ::
     |  nmlat = 97,         ! number of magnetic latitudes
     |  nmlon = 80,         ! number of magnetic longitudes
     |  nmlonp1=nmlon+1,
     |  nmlonp2=nmlon+2,
     |  nmlatp1=nmlat+1,
     |  nmlath=(nmlat+1)/2, ! index to magnetic equator
     |  nmlev   = nlev,     ! number of midpoint levels of mag fields (1,nlev)
     |  nimlev  = nmlev,    ! number of interface levels of mag fields (1,nlev)
     |  nmlevp1 = nlevp1,   ! 
     |  nimlevp1= nimlev+1  ! 
!
! kbotdyn is zp k-index at which upward dynamo integrals begin (dynamo.F):
! For zp -8 at dzp=0.5,  kbotdyn = 5 
! For zp -8 at dzp=0.25, kbotdyn = 9 
!
! -10 to +7 by 0.50 ("normal" vertical resolution)
#if (NLEV==34)   
!     integer,parameter :: kbotdyn = 1
!
! -10 to +7 by 0.25 ("double" vertical resolution)
#elif (NLEV==68) 
!     integer,parameter :: kbotdyn = 1
#endif
!
! For modules hist_mod and input_mod:
      character(len=16),parameter :: tgcm_version = 'vtgcm_trunk     '
!
! tgcm_name is either "tiegcm" or "time-gcm"
      character(len=8),parameter :: tgcm_name    = 'vtgcm'
!
      integer,parameter ::
     |  mxday=366,         ! max allowed calendar day
     |  mxhvols = 100,     ! max number of output history files
     |  mxseries = 10,     ! max number of time series for primary histories
     |  mxseries_sech = 50,! max number of secondary history time series
     |  mxfsech = 100,     ! max number of fields on secondary histories
     |  mxind_time = 100,  ! max number of time-dependent solar index points
!
! For now (3/26/13), mxnlev_oxvgcm must be the same as nlev in the oxvgcm_ncfile
     |  mxnlev_oxvgcm = 1  ! max number of levels to nudge with oxvgcm
!
      end module params_module
