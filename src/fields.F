!
      module fields_module
!
! There are allocatable arrays and an array of field structures for 
!   3-d and 4-d fields. Subdomains are allocated to the allocatable 
!   arrays, and structure data pointers point to the allocatable arrays,
!   e.g. f4d(i_tn)%data(:,:,:,:) => tn(:,:,:,:), where the 4 dimensions
!   are (pressure,longitude,latitude,2). The last dimension is for
!   previous and current time steps. 3-d fields are the same as 4-d
!   without the final time dimension. 
!   -- latest revision: (12/04/06) swb: add 3-D fields (difk,dift,xmue)
!   -- latest revision: (03/31/08) swb: use 4-D field (n4s_nm) for n4s
!                                       replaces psb_nm, i_psb_nm
!   -- latest revision: (05/12/08) swb: use 4-D field (no, no_nm) for no
!                                       replaces psc_nm, i_psc_nm
!                                       replaces i_eno
!   -- current revision  : (09/24/13) bf : use 4-D field (co_nm) for co
!                                       replaces psa_nm, i_psa_nm
! btf 6/13/13: 
!   - Add 4d prognostics so, so_nm, so2, so2_nm, incrementing nf4d from 25 to 29
!   - Add 3d fields vso,vso2,vso3,scso,scso2,scso3 incrementing nf3d from 34 to 40
! btf 9/24/13:
!   - Replace psa_nm with co_nm
!
      use params_module,only: nlonp1,nlonp4,nlat,nlatp1,nlevp1,
     |  mxfsech,spval,nmlat,nmlonp1
      implicit none
      save
      integer,parameter :: 
     |  nf4d          = 29, ! number of 4-d fields
     |  nf4d_hist     = 28, ! number of 4-d fields on primary histories
     |  nf3d          = 40  ! number of 3-d fields
      real,parameter :: 
     |  field_initval = 0.  ! initialization value for fields data
      integer,parameter :: 
     |  longname_len  = 80, ! length of field long name 
     |  shortname_len = 16, ! length of field short name 
     |  units_len     = 16  ! length of field units attribute
      integer :: levd0,levd1, lond0,lond1, latd0,latd1 ! subdomain dimensions
!
! 4-d fields structure type:
! Data pointers will point to allocatable target arrays.
!
      type fields_4d
        character(len=longname_len)  :: long_name
        character(len=shortname_len) :: short_name
        character(len=units_len)     :: units
        character(len=16)            :: vcoord ! midpoints or interfaces
        logical                      :: 
     |    magnetic,    ! true if field is on magnetic 3d grid
     |    prognostic,  ! true if prognostic (diagnostic if false)
     |    mpi          ! flag used by some mpi routines
        real                         :: polesign
        real,pointer                 :: data(:,:,:,:)
      end type fields_4d
!
! 21 4-d fields for vtgcm:
! tn,un,vn,o1,co,n4s,n2d,no,eno,o2,n2,w,z,
! tn_nm,un_nm,vn_nm,o2_nm,o1_nm,co_nm,n4s_nm,no_nm,psh_nm
!
! Original vtgcm ndex.h:
!                                    o   co
!     COMMON/NDEX/NDEXA(1),NT,NU,NV,NPS,NPS2,NPSA,NPSB,NPSC,NPSH,NW,NZ,
!    .NTNM,NUNM,NVNM,NPSNM,NPS2NM,NPSANM,NPSBNM,NPSCNM,NPSHNM,NMS,
!    .NVC,NKLDU,NKLDV,NKLDT,NKLDPS,NKLDP2,NKLDPA,NKLDPB,NKLDPC,NKLDPH,
!    .NDEXB(1),NRH,NFLH,NFPH,NQDH,NKMH,NPSDH,NPSDH2,NPSADH,NPSBDH,
!    .NPSCDH,NPSHDH,NLXX,NLYY,NLXY,NKT,NKM,NCP,NQI,NRJI,
!    .NDEXC(1),NQ,NUI,NVI,NRJ,NQUV,NDEXD
!
      type (fields_4d) :: f4d(nf4d)
!
! Indices, to f4d fields, e.g. f4d(i_tn)%data is neutral temperature.
      integer ::
     |  i_tn    ,i_un     ,i_vn     ,i_o1     ,i_co     ,
     |  i_n4s   ,i_n2d    ,i_no     ,i_o2     ,i_n2     ,
     |  i_w     ,i_z      ,i_tn_nm  ,i_un_nm  ,i_vn_nm  ,
     |  i_o2_nm ,i_o1_nm  ,i_co_nm  ,i_n4s_nm ,i_no_nm  ,
     |  i_psh_nm,i_o2p    ,i_op     ,i_ne     ,i_so     ,
     |  i_so_nm  ,i_so2    ,i_so2_nm ,i_vc
!
! Allocatable target arrays for 4-d fields:
! Allocation will be:   tn(levd0:levd1,lond0 :lond1 ,latd0 :latd1 ,2)
! which is the same as: tn(1    :nlev ,lon0-2:lon1+2,lat0-2:lat1+2,2)
!
      real,dimension(:,:,:,:),target,allocatable :: ! (k,i,j,2)
     |  tn      ,un      ,vn      ,o1      ,co      ,
     |  n4s     ,n2d     ,no      ,o2      ,n2      ,
     |  w       ,z       ,tn_nm   ,un_nm   ,vn_nm   ,
     |  o2_nm   ,o1_nm   ,co_nm   ,n4s_nm  ,no_nm   ,
     |  psh_nm  ,o2p     ,op      ,ne      ,so      ,
     |  so_nm   ,so2     ,so2_nm  ,vc
!
! 3-d fields structure type:
! Data pointers will point to allocatable target arrays.
!
      type fields_3d
        character(len=longname_len)  :: long_name
        character(len=shortname_len) :: short_name
        character(len=units_len)     :: units
!       character(len=16)            :: vcoord ! midpoints or interfaces
        character(len=19)            :: vcoord ! midpoints or interfaces
        character(len=8)             :: dimnames(3) ! e.g., lat,lon,lev, or lat,lon,' '
        logical                      :: 
     |    magnetic,    ! true if field is on magnetic 3d grid 
     |    prognostic,  ! true if prognostic (diagnostic if false)
     |    mpi,         ! flag used by some mpi routines
     |    task0_only   ! if true, field was defined at root task only
        real                         :: polesign
        real,pointer                 :: data(:,:,:)
      end type fields_3d
      type (fields_3d) :: f3d(nf3d)
!
! 3-d fields and indices, e.g., reference f3d(kldt)%data
      integer :: 
     |  i_kldt   ,i_kldu   ,i_kldv    ,i_kldo2 ,i_kldo1 ,
     |  i_kldco  ,i_barm   ,i_cp      ,i_kt    ,
     |  i_km     ,i_hdt    ,i_hdu     ,i_hdv   ,i_hdo1  ,
     |  i_hdco   ,i_vo1    ,i_vco     ,i_vco2  ,i_vn2   ,
     |  i_vo2    ,i_sco1   ,i_scco    ,i_scco2 ,i_scn2  ,
     |  i_sco2   ,i_co2    ,i_n2p     ,i_nop   ,i_co2p  ,
     |  i_te     ,i_ti     ,i_difk    ,i_dift  ,i_xmue  ,
     |  i_vso    ,i_vso2   ,i_vso3    ,i_scso  ,i_scso2 ,
     |  i_scso3
!
! Allocatable target arrays 3-d:
! Allocation will be:   tn(levd0:levd1,lond0 :lond1 ,latd0 :latd1 )
! which is the same as: tn(1    :nlev ,lon0-2:lon1+2,lat0-2:lat1+2)
!
      real,dimension(:,:,:),target,allocatable :: ! (k,i,j)
     |  kldt     ,kldu     ,kldv     ,kldo2    ,kldo1   ,
     |  kldco    ,barm     ,cp       ,kt       ,
     |  km       ,hdt      ,hdu      ,hdv      ,hdo1    ,
     |  hdco     ,vo1      ,vco      ,vco2     ,vn2     ,
     |  vo2      ,sco1     ,scco     ,scco2    ,scn2    ,
     |  sco2     ,co2      ,n2p      ,nop      ,co2p    ,
     |  te       ,ti       ,difk     ,dift     ,xmue    ,
     |  vso      ,vso2     ,vso3     ,scso     ,scso2   ,
     |  scso3
!
! 2-d field type (used for 2d secondary history fields):
      type fields_2d
        character(len=longname_len)  :: long_name
        character(len=shortname_len) :: short_name
        character(len=units_len)     :: units
        logical                      :: 
     |    magnetic,    ! true if field is on magnetic 3d grid 
     |    prognostic,  ! true if prognostic (diagnostic if false)
     |    mpi,         ! flag used by some mpi routines
     |    task0_only   ! if true, field was defined at root task only
        real                         :: polesign
        real,pointer                 :: data(:,:) ! (k,i)
      end type fields_2d
!
! Electric potential on geographic and magnetic grids:
! (full domains until dynamo is parallelized)
      real ::
     |  dynpot(nlonp1,0:nlatp1,nlevp1),  ! 3d electric potential geographic
     |  phim3d(nmlonp1,nmlat,nlevp1),    ! 3d electric potential magnetic
     |  emphi3d(nmlonp1,nmlat,nlevp1),   ! 3d eastward electric field magnetic
     |  emlam3d(nmlonp1,nmlat,nlevp1),   ! 3d equatorw. electric field magnetic
     |  emz3d(nmlonp1,nmlat,nlevp1)      ! 3d upward (?) electric field magnetic
!
! Secondary history field structures.
! These are initialized in sub init_fsech, and data is set in addfld.F.
      type fields_sech
        character(len=longname_len)  :: long_name
        character(len=shortname_len) :: short_name
        character(len=units_len)     :: units
        logical                      :: task0_only
        logical                      :: prognostic
        logical                      :: mag,geo
        integer                      :: ndims        ! 2d or 3d
        character(len=8)             :: dimnames(3)  ! dims on history
        integer                      :: dimsizes(3)  ! dim sizes on history
        real,pointer                 :: data(:,:,:)  ! allocated by addfld
      end type fields_sech
      type(fields_sech) :: fsechist(mxfsech)
!
! Time indices for rightmost dimension of 4d data, itp for current timestep, 
! itc for next timestep. Fields at the previous timestep (time n-1) are
! saved at both itp and itc (e.g., tn_nm, un_nm, etc).
!
      integer :: itc,itp 
!
! If fakeflds is true, use fake dimensions for fields
! (for dry runs, testing, etc)
!
      logical,parameter :: fakeflds = .false.
!
      real,dimension(nlevp1,nlonp4+1,-2:nlat) ::
     |  fnrh,  ! eddy viscosity
     |  fkmh   ! M/T
!
! Lower boundary conditions (bottom interface level) for t,u,v from 
! source history, and at current model time:
      real,allocatable,dimension(:,:) ::
     |  tlbc,     ulbc,     vlbc,     ! subdomains (lond0:lond1,latd0:latd1)
     |  tlbc_nm,  ulbc_nm,  vlbc_nm,  ! subdomains (lond0:lond1,latd0:latd1)
     |  tlbc_glb, ulbc_glb, vlbc_glb, ! global (nlonp4,nlat)
     |  tlbc_nm_glb, ulbc_nm_glb, vlbc_nm_glb  ! global (nlonp4,nlat)
!
! Full 3d grid with all primary history fields for writing to netcdf
! history files. This will be allocated only on the root task 
! (see allocdata):
      real,allocatable :: foutput(:,:,:,:) ! (nlevp1,nlonp4,nlat,nf4d_hist)
      real,allocatable :: fzg(:,:,:)       ! (nlevp1,nlonp4,nlat)
!
! Parameters declared to maintain skeletal code for currently
! unused modules (e.g., dynamo)
      integer :: i_poten=0
      contains
!-----------------------------------------------------------------------
!
      subroutine init_4d(lon0,lon1,lat0,lat1,mytid,iprint)
!
! Set names, units, indices and pointers for f4d and f3d field structures,
! and allocate 3d and 4d field arrays. Also make other data allocations.
!
      implicit none
!
! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1,mytid,iprint
!
! Local:
      integer :: n,istat,i
!
! Fields are allocated at full task subdomain, including ghost cells:
      levd0 = 1      ; levd1 = nlevp1
      lond0 = lon0-2 ; lond1 = lon1+2
      latd0 = lat0-2 ; latd1 = lat1+2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      n = 1
!
! 4-d fields (long and short names, units, index):
! (pointer definition must follow allocate statement)
! Note it is difficult to make a subroutine to do the allocations because
!   allocatable arrays cannot be dummy arguments. (Once they are allocated,
!   they can be passed to subroutines and declared as real, but subroutine 
!   dummy arguments cannot be declared as allocatable arrays.)
!
! vtgcm 4d:
!     real,dimension(:,:,:,:),target,allocatable :: ! (k,i,j,2)
!    |  tn      ,un      ,vn      ,o1      ,co      ,
!    |  n4s     ,n2d     ,eno     ,o2      ,n2      ,
!    |  w       ,z       ,tn_nm   ,un_nm   ,vn_nm   ,
!    |  o2_nm   ,o1_nm   ,co_nm   ,n4s_nm  ,psc_nm  ,
!    |  psh_nm  ,op2     ,op      ,ne      ,vc

      allocate(tn(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NEUTRAL TEMPERATURE')")
      f4d(n)%short_name = "TN"
      f4d(n)%units = "K"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => tn
      i_tn = n ; n = n+1
!
      allocate(un(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NEUTRAL ZONAL WIND (+EAST)')")
      f4d(n)%short_name = "UN"
      f4d(n)%units = "cm/s"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => un
      i_un = n ; n = n+1
!
      allocate(vn(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NEUTRAL MERIDIONAL WIND (+NORTH)')")
      f4d(n)%short_name = "VN"
      f4d(n)%units = "cm/s"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => vn
      i_vn = n ; n = n+1
!
      allocate(o1(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('ATOMIC OXYGEN')")
      f4d(n)%short_name = ("O1")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => o1
      i_o1 = n ; n = n+1
!
      allocate(co(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('CARBON MONOXIDE')")
      f4d(n)%short_name = ("CO")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => co
      i_co = n ; n = n+1
!
      allocate(n4s(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('N4S')")
      f4d(n)%short_name = ("N4S")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => n4s
      i_n4s = n ; n = n+1
!
      allocate(n2d(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('N2D')")
      f4d(n)%short_name = ("N2D")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => n2d
      i_n2d = n ; n = n+1
!
      allocate(no(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NO')")
      f4d(n)%short_name = ("NO")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => no
      i_no = n ; n = n+1
!
      allocate(o2(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('MOLECULAR OXYGEN')")
      f4d(n)%short_name = ("O2")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => o2
      i_o2 = n ; n = n+1
!
      allocate(n2(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('MOLECULAR NITROGEN')")
      f4d(n)%short_name = ("N2")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => n2
      i_n2 = n ; n = n+1
!
      allocate(w(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('VERTICAL MOTION')")
      f4d(n)%short_name = ("OMEGA")
      f4d(n)%units = "s-1"
      f4d(n)%vcoord = "interfaces"
      f4d(n)%data => w
      i_w = n ; n = n+1
!
      allocate(z(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('GEOPOTENTIAL HEIGHT')")
      f4d(n)%short_name = ("Z")
      f4d(n)%units = "cm"
      f4d(n)%vcoord = "interfaces"
      f4d(n)%data => z
      i_z = n ; n = n+1
!
      allocate(tn_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NEUTRAL TEMPERATURE (TIME N-1)')")
      f4d(n)%short_name = ("TN_NM")
      f4d(n)%units = "K"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => tn_nm
      i_tn_nm = n ; n = n+1
!
      allocate(un_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NEUTRAL ZONAL WIND (TIME N-1)')")
      f4d(n)%short_name = ("UN_NM")
      f4d(n)%units = "cm/s"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => un_nm
      i_un_nm = n ; n = n+1
!
      allocate(vn_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NEUTRAL MERIDIONAL WIND (TIME N-1)')")
      f4d(n)%short_name = ("VN_NM")
      f4d(n)%units = "cm/s"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => vn_nm
      i_vn_nm = n ; n = n+1
!
      allocate(o2_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('MOLECULAR OXYGEN (TIME N-1)')")
      f4d(n)%short_name = ("O2_NM")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => o2_nm
      i_o2_nm = n ; n = n+1
!
      allocate(o1_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('ATOMIC OXYGEN (TIME N-1)')")
      f4d(n)%short_name = ("O1_NM")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => o1_nm
      i_o1_nm = n ; n = n+1
!
      allocate(co_nm (levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('CARBON MONOXIDE (TIME N-1)')")
      f4d(n)%short_name = ("CO_NM")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => co_nm 
      i_co_nm  = n ; n = n+1
!
      allocate(n4s_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('N4S (TIME N-1)')")
      f4d(n)%short_name = ("N4S_NM")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => n4s_nm
      i_n4s_nm = n ; n = n+1
!
      allocate(no_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NO (TIME N-1)')")
      f4d(n)%short_name = ("NO_NM")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => no_nm
      i_no_nm = n ; n = n+1
!
      allocate(psh_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('PSH (TIME N-1)')")
      f4d(n)%short_name = ("PSH_NM")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => psh_nm
      i_psh_nm = n ; n = n+1
!
      allocate(o2p(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('O2+')")
      f4d(n)%short_name = ("O2P")
      f4d(n)%units = "cm-3"
      f4d(n)%vcoord = "interfaces"
      f4d(n)%data => o2p
      i_o2p = n ; n = n+1
!
      allocate(op(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('O+')")
      f4d(n)%short_name = ("OP")
      f4d(n)%units = "cm-3"
      f4d(n)%vcoord = "interfaces"
      f4d(n)%data => op
      i_op = n ; n = n+1
!
      allocate(ne(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('ELECTRON DENSITY')")
      f4d(n)%short_name = ("NE")
      f4d(n)%units = "cm-3"
      f4d(n)%vcoord = "interfaces"
      f4d(n)%data => ne
      i_ne = n ; n = n+1
!
      allocate(so(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('SULPHUR MONOXIDE')")
      f4d(n)%short_name = ("SO")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => so
      i_so = n ; n = n+1
!
      allocate(so_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('SULPHUR MONOXIDE (TIME N-1)')")
      f4d(n)%short_name = ("SO_NM")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => so_nm
      i_so_nm = n ; n = n+1
!
      allocate(so2(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('SULPHUR DIOXIDE')")
      f4d(n)%short_name = ("SO2")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => so2
      i_so2 = n ; n = n+1
!
      allocate(so2_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('SULPHUR DIOXIDE (TIME N-1)')")
      f4d(n)%short_name = ("SO2_NM")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => so2_nm
      i_so2_nm = n ; n = n+1
!
      allocate(vc(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('COS(PHI)*VN')")
      f4d(n)%short_name = ("VC")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => vc
      i_vc = n ; n = n+1
!
      if (n-1 /= nf4d) then
        write(6,"('>>> fields: n=',i3,' nf3d=',i3)") n,nf4d
        call shutdown('nf4d')
      endif
!
! itp,itc are indices to rightmost dimension of field data, for
! previous and current time steps.
      itp = 1
      itc = 2
!
! Set polesign for crossing the poles (see mk_polelats and mp_bndlats)
!
! 11/05 btf: Getting the following error here, attempting to initialize
!            f4d(:)%polesign, from hao (callisto) pgf90 5.1-6:
!
! Lowering Error: symbol data$sd is a member reference
! Lowering Error: symbol data$sd is a member reference
! PGF90-F-0000-Internal compiler error. Errors in Lowering      21 (/home/foster/tgcm/tiegcm-hist-lbc/src/fields.F: 437)
! PGF90/any Linux/x86 5.1-6: compilation aborted
!
! However, this error does not occur if fields.F is local!
!
!     f4d(:)%polesign = 1. ! 11/05 btf: pgf90 at hao does not like this
!     do i=1,size(f4d)
      do i=1,nf4d
        f4d(i)%polesign = 1.
      enddo
      f4d(i_un)%polesign = -1.
      f4d(i_un_nm)%polesign = -1.
      f4d(i_vn)%polesign = -1.
      f4d(i_vn_nm)%polesign = -1.
      f4d(i_n2d)%polesign = 0.
!     f4d(i_ne )%polesign = 0.
!     f4d(i_o2p)%polesign = 0.
!
! All f4d fields are on geographic grid:
      f4d%magnetic = .false.
!
! f4d fields are "prognostic":
      f4d%prognostic = .true.
!
! Data init:
      do i=1,nf4d
        f4d(i)%data = field_initval
      enddo
!
! Report to stdout:
      if (iprint > 0) then
        write(6,"(/,'There are ',i3,' 4-d fields:')") nf4d
        do n=1,nf4d
          call print_f4d(f4d(n),n)
        enddo
        write(6,"('Initialized f4d%data to ',e12.4)") field_initval
      endif
      end subroutine init_4d
!-----------------------------------------------------------------------
      subroutine init_3d(lon0,lon1,lat0,lat1,mytid,iprint)
!
! Set names, units, indices and pointers for f4d and f3d field structures,
! and allocate 3d and 4d field arrays. Also make other data allocations.
!
      implicit none
!
! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1,mytid,iprint
!
! Local:
      integer :: i,n,istat
!
      n = 1
!
! 3-d fields (long and short names, units, index):
!
      allocate(kldt(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('TN HORIZONTAL DIFFUSION COEFF')")
      f3d(n)%short_name = "KLDT"
      f3d(n)%units = " "
      f3d(n)%data => kldt
      i_kldt = n ; n = n+1
!
      allocate(kldu(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('UN HORIZONTAL DIFFUSION COEFF')")
      f3d(n)%short_name = "KLDU"
      f3d(n)%units = " "
      f3d(n)%data => kldu
      i_kldu = n ; n = n+1
!
      allocate(kldv(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('VN HORIZONTAL DIFFUSION COEFF')")
      f3d(n)%short_name = "KLDV"
      f3d(n)%units = " "
      f3d(n)%data => kldv
      i_kldv = n ; n = n+1
!
      allocate(kldo2(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('O2 HORIZONTAL DIFFUSION COEFF')")
      f3d(n)%short_name = "KLDO2"
      f3d(n)%units = " "
      f3d(n)%data => kldo2
      i_kldo2 = n ; n = n+1
!
      allocate(kldo1(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('O1 HORIZONTAL DIFFUSION COEFF')")
      f3d(n)%short_name = "KLDO1"
      f3d(n)%units = " "
      f3d(n)%data => kldo1
      i_kldo1 = n ; n = n+1
!
      allocate(kldco(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('CO HORIZONTAL DIFFUSION COEFF')")
      f3d(n)%short_name = "KLDCO"
      f3d(n)%units = " "
      f3d(n)%data => kldco
      i_kldco = n ; n = n+1
!
      allocate(barm(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('MEAN MOLECULAR WEIGHT (MBAR)')")
      f3d(n)%short_name = "BARM"
      f3d(n)%units = " "
      f3d(n)%data => barm
      i_barm = n ; n = n+1
!
      allocate(cp(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('SPECIFIC HEAT (CP)')")
      f3d(n)%short_name = "CP"
      f3d(n)%units = " "
      f3d(n)%data => cp
      i_cp = n ; n = n+1
!
      allocate(kt(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('MOLECULAR THERMAL CONDUCTIVITY (KT)')")
      f3d(n)%short_name = "KT"
      f3d(n)%units = " "
      f3d(n)%data => kt
      i_kt = n ; n = n+1
!
      allocate(km(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('MOLECULAR DIFFUSION (KM)')")
      f3d(n)%short_name = "KM"
      f3d(n)%units = " "
      f3d(n)%data => km
      i_km = n ; n = n+1
!
      allocate(hdt(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('TN HORIZONTAL DIFFUSION')")
      f3d(n)%short_name = "HDT"
      f3d(n)%units = " "
      f3d(n)%data => hdt
      i_hdt = n ; n = n+1
!
      allocate(hdu(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('UN HORIZONTAL DIFFUSION')")
      f3d(n)%short_name = "HDU"
      f3d(n)%units = " "
      f3d(n)%data => hdu
      i_hdu = n ; n = n+1
!
      allocate(hdv(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('VN HORIZONTAL DIFFUSION')")
      f3d(n)%short_name = "HDV"
      f3d(n)%units = " "
      f3d(n)%data => hdv
      i_hdv = n ; n = n+1
!
      allocate(hdo1(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('O1 HORIZONTAL DIFFUSION')")
      f3d(n)%short_name = "HDO1"
      f3d(n)%units = " "
      f3d(n)%data => hdo1
      i_hdo1 = n ; n = n+1
!
      allocate(hdco(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('CO HORIZONTAL DIFFUSION')")
      f3d(n)%short_name = "HDCO"
      f3d(n)%units = " "
      f3d(n)%data => hdco
      i_hdco = n ; n = n+1
!
      allocate(vo1(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('O1 LINE INTEGRAL')")
      f3d(n)%short_name = "VO1"
      f3d(n)%units = " "
      f3d(n)%data => vo1
      i_vo1 = n ; n = n+1
!
      allocate(vco(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('CO LINE INTEGRAL')")
      f3d(n)%short_name = "VCO"
      f3d(n)%units = " "
      f3d(n)%data => vco
      i_vco = n ; n = n+1
!
      allocate(vco2(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('CO2 LINE INTEGRAL')")
      f3d(n)%short_name = "VCO2"
      f3d(n)%units = " "
      f3d(n)%data => vco2
      i_vco2 = n ; n = n+1
!
      allocate(vn2(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('N2 LINE INTEGRAL')")
      f3d(n)%short_name = "VN2"
      f3d(n)%units = " "
      f3d(n)%data => vn2
      i_vn2 = n ; n = n+1
!
      allocate(vo2(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('O2 LINE INTEGRAL')")
      f3d(n)%short_name = "VO2"
      f3d(n)%units = " "
      f3d(n)%data => vo2
      i_vo2 = n ; n = n+1
!
      allocate(sco1(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('O1 SLANT COLUMN DENSITY')")
      f3d(n)%short_name = "VO1"
      f3d(n)%units = " "
      f3d(n)%data => sco1
      i_sco1 = n ; n = n+1
!
      allocate(scco(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('CO SLANT COLUMN DENSITY')")
      f3d(n)%short_name = "VCO"
      f3d(n)%units = " "
      f3d(n)%data => scco
      i_scco = n ; n = n+1
!
      allocate(scco2(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('CO2 SLANT COLUMN DENSITY')")
      f3d(n)%short_name = "SCCO2"
      f3d(n)%units = " "
      f3d(n)%data => scco2
      i_scco2 = n ; n = n+1
!
      allocate(scn2(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('N2 SLANT COLUMN DENSITY')")
      f3d(n)%short_name = "SCN2"
      f3d(n)%units = " "
      f3d(n)%data => scn2
      i_scn2 = n ; n = n+1
!
      allocate(sco2(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('O2 SLANT COLUMN DENSITY')")
      f3d(n)%short_name = "SCO2"
      f3d(n)%units = " "
      f3d(n)%data => sco2
      i_sco2 = n ; n = n+1
!
      allocate(co2(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('CARBON DIOXIDE (CO2)')") 
      f3d(n)%short_name = "CO2"
      f3d(n)%units = " "
      f3d(n)%data => co2
      i_co2 = n ; n = n+1
!
      allocate(n2p(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('N2+')") 
      f3d(n)%short_name = "N2P"
      f3d(n)%units = " "
      f3d(n)%data => n2p
      i_n2p = n ; n = n+1
!
      allocate(nop(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('NO+')") 
      f3d(n)%short_name = "NOP"
      f3d(n)%units = " "
      f3d(n)%data => nop
      i_nop = n ; n = n+1
!
      allocate(co2p(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('CO2+')") 
      f3d(n)%short_name = "CO2P"
      f3d(n)%units = " "
      f3d(n)%data => co2p
      i_co2p = n ; n = n+1
!
      allocate(te(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('ELECTRON TEMPERATURE')") 
      f3d(n)%short_name = "TE"
      f3d(n)%units = " "
      f3d(n)%data => te
      i_te = n ; n = n+1
!
      allocate(ti(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('ION TEMPERATURE')") 
      f3d(n)%short_name = "TI"
      f3d(n)%units = " "
      f3d(n)%data => ti
      i_ti = n ; n = n+1
!
      allocate(difk(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('EDDY DIFFUSION COEFFICIENT')") 
      f3d(n)%short_name = "DIFK"
      f3d(n)%units = "(1/sec)"
      f3d(n)%data => difk
      i_difk = n ; n = n+1
!
      allocate(dift(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('EDDY CONDUCTIVITY COEFFICIENT')") 
      f3d(n)%short_name = "DIFT"
      f3d(n)%units = "(1/sec)"
      f3d(n)%data => dift
      i_dift = n ; n = n+1
!
      allocate(xmue(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('EDDY VISCOSITY COEFFICIENT')") 
      f3d(n)%short_name = "XMUE"
      f3d(n)%units = "(1/sec)"
      f3d(n)%data => xmue
      i_xmue = n ; n = n+1
!
      allocate(vso(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('SO LINE INTEGRAL ')")
      f3d(n)%short_name = "VSO"
      f3d(n)%units = " "
      f3d(n)%data => vso
      i_vso = n ; n = n+1
!
      allocate(vso2(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('SO2 LINE INTEGRAL ')")
      f3d(n)%short_name = "VSO2"
      f3d(n)%units = " "
      f3d(n)%data => vso2
      i_vso2 = n ; n = n+1
!
      allocate(vso3(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('SO3 LINE INTEGRAL ')")
      f3d(n)%short_name = "VSO3"
      f3d(n)%units = " "
      f3d(n)%data => vso3
      i_vso3 = n ; n = n+1
!
      allocate(scso(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('SO SLANT COLUMN DENSITY ')")
      f3d(n)%short_name = "SCSO"
      f3d(n)%units = " "
      f3d(n)%data => scso
      i_scso = n ; n = n+1
!
      allocate(scso2(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('SO2 SLANT COLUMN DENSITY ')")
      f3d(n)%short_name = "SCSO2"
      f3d(n)%units = " "
      f3d(n)%data => scso2
      i_scso2 = n ; n = n+1
!
      allocate(scso3(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('SO3 SLANT COLUMN DENSITY ')")
      f3d(n)%short_name = "SCSO3"
      f3d(n)%units = " "
      f3d(n)%data => scso3
      i_scso3 = n ; n = n+1
!
      if (n-1 /= nf3d) then
        write(6,"('>>> fields: n=',i3,' nf3d=',i3)") n,nf3d
        call shutdown('nf3d')
      endif
!
! f3d fields are diagnostics:
      f3d%prognostic = .false.
      f3d%magnetic = .false.
      f3d%task0_only = .false.
      f3d%vcoord = 'midpoints'
      do n=1,nf3d
        if (trim(f3d(n)%short_name)=='ZG') f3d(n)%vcoord='interfaces' 
      enddo
!
! Data init:
      do i=1,nf3d
        f3d(i)%data = field_initval
      enddo
!
! Report to stdout:
      if (iprint > 0) then
        write(6,"(/,'There are ',i3,' 3-d fields:')") nf3d
        do n=1,nf3d
          call print_f3d(f3d(n),n)
        enddo
      endif
      end subroutine init_3d
!-----------------------------------------------------------------------
      subroutine init_fsech
!
! Initialize secondary history field structures. This does not include
! data, which is allocated and defined by sub addfld (addfld.F).
!
      use input_module,only: secflds
      use hist_module,only: nfsech
      implicit none
!
! Local:
      integer :: i,iprog,ier
!
! External:
      integer,external :: strloc
!
! Secondary history fields:
      fsechist(:)%long_name = ' '
      fsechist(:)%short_name= ' '
      fsechist(:)%units     = ' '
      fsechist(:)%prognostic  = .false.
      fsechist(:)%task0_only  = .false.
      fsechist(:)%mag = .false.
      fsechist(:)%geo = .false.

      do i=1,nfsech
        if (len_trim(secflds(i)) > 0) then
          iprog = strloc(f4d%short_name,nf4d,secflds(i))
!
! Is a prognostic: define fsech(i) from f4d(iprog):
!
          if (iprog > 0) then
            fsechist(i)%prognostic = .true.
            fsechist(i)%short_name = f4d(iprog)%short_name
            fsechist(i)%long_name  = f4d(iprog)%long_name
            fsechist(i)%units      = f4d(iprog)%units
!
! Is a diagnostic: define fsech(i)%name from input field name.
! Set units and long_name blank (will be optionally defined in
! user called sub addfsech)
!
          else ! is diagnostic
            fsechist(i)%long_name = secflds(i)
            fsechist(i)%short_name = secflds(i)
            fsechist(i)%units = ' '
            fsechist(i)%long_name = ' '
          endif
        endif
      enddo ! i=1,nfsech
      end subroutine init_fsech
!-----------------------------------------------------------------------
      subroutine init_lbc(lon0,lon1,lat0,lat1)
      implicit none
!
! Arg:
      integer,intent(in) :: lon0,lon1,lat0,lat1
!
! Local:
      integer :: istat
!
! Allocate t,u,v lbc (t,u,v at bottom interface level)
! These will be read from source history, then set in dt.F and duv.F.
!
! Subdomains:
      allocate(tlbc(lond0:lond1,latd0:latd1),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of tlbc')
      allocate(ulbc(lond0:lond1,latd0:latd1),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of ulbc')
      allocate(vlbc(lond0:lond1,latd0:latd1),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of vlbc')
      write(6,"('init_lbc: allocated subdomains tlbc, ulbc, vlbc')")
      
      allocate(tlbc_nm(lond0:lond1,latd0:latd1),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of tlbc_nm')
      allocate(ulbc_nm(lond0:lond1,latd0:latd1),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of ulbc_nm')
      allocate(vlbc_nm(lond0:lond1,latd0:latd1),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of vlbc_nm')
      write(6,"('init_lbc_nm: allocated subdomains tlbc_nm, ulbc_nm,',
     |  ' vlbc_nm')")
!
! Global domain (subdomains will be gathered into these arrays
! for writing to the history, see sub mp_gather2root_lbc in mpi.F):
      allocate(tlbc_glb(nlonp4,nlat),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of tlbc_glb')
      allocate(ulbc_glb(nlonp4,nlat),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of ulbc_glb')
      allocate(vlbc_glb(nlonp4,nlat),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of vlbc_glb')
      write(6,"('init_lbc: allocated globals tlbc_glb, ulbc_glb,',
     |  ' vlbc_glb: nlonp4=',i3,' nlat=',i3)") nlonp4,nlat
     
      allocate(tlbc_nm_glb(nlonp4,nlat),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of tlbc_nm_glb')
      allocate(ulbc_nm_glb(nlonp4,nlat),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of ulbc_nm_glb')
      allocate(vlbc_nm_glb(nlonp4,nlat),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of vlbc_nm_glb')
      write(6,"('init_lbc: allocated globals tlbc_nm_glb, ulbc_nm_glb,',
     |  ' vlbc_nm_glb: nlonp4=',i3,' nlat=',i3)") nlonp4,nlat
!
      end subroutine init_lbc
!-----------------------------------------------------------------------
      subroutine print_f4d(f,n)
      implicit none
      type(fields_4d),intent(in) :: f
      integer,intent(in) :: n
!
      write(6,"('Field ',i3,': ',a,' Short name: ',a,' Units: ',a)") 
     |  n,f%long_name(1:40),f%short_name(1:8),f%units(1:8)
      end subroutine print_f4d
!-----------------------------------------------------------------------
      subroutine print_f3d(f,n)
      implicit none
      type(fields_3d),intent(in) :: f
      integer,intent(in) :: n
!
      write(6,"('Field ',i3,': ',a,' Short name: ',a,' Units: ',a)") 
     |  n,f%long_name(1:40),f%short_name(1:8),f%units(1:8)
      end subroutine print_f3d
!-----------------------------------------------------------------------
      subroutine print_f2d(f,n)
      implicit none
      type(fields_2d),intent(in) :: f
      integer,intent(in) :: n
!
      write(6,"('Field ',i3,': ',a,' Short name: ',a,' Units: ',a)") 
     |  n,f%long_name(1:40),f%short_name(1:8),f%units(1:8)
      end subroutine print_f2d
!-----------------------------------------------------------------------
      end module fields_module
