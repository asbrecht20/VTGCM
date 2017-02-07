      module oxvgcm
!
! Read lower boundary conditions from Oxford Venus GCM (Chris Lee), 
! and use them to nudge vtgcm t,u,v from its bottom boundary upward.
! (lower boundary of Z is set from the data, but not nudged upward)
!
! As of Nov 20, 2012, we have a set of lower boundary data for t,uj,v,z
! at 10 zp levels (-16 to -11.5 by 0.5) at a single time.
!
! If oxvgcm_zm==1, then zonally average before nudging.
!
      use params_module,only: nlon,nlonp4,nlat,nlevp1,dlev,zibot,
     | mxnlev_oxvgcm
      use input_module ,only: 
     |  oxvgcm_ncfile,        ! Oxford VGCM data file
     |  oxvgcm_zm             ! Flag to nudge with zonal means
      use cons_module  ,only: pi
      use init_module  ,only: istep
      use netcdf
      implicit none
!#include <netcdf.inc>
      save
      integer,parameter ::  
     |  mxtime=20  ! max number of times allowed on file

      integer :: ncid=0   ! netcdf file id
      integer :: ntime=0  ! number of times on the file
!
! Number of levels to nudge upward from bottom boundary (zp -16):
!     integer,parameter :: oxvgcm_nlev = 10  ! nudge from -16 to zp -11.5
!     integer,parameter :: oxvgcm_nlev = 3  ! nudge from -16 to zp -15.0
      integer,parameter :: oxvgcm_nlev = 1  ! nudge from -16  (only)
!
! Lower boundaries from oxvgcm_ncfile:
      real,dimension(nlonp4,nlat,mxnlev_oxvgcm) :: ! lbc for t,u,v,z
     |  t_oxvgcm,u_oxvgcm,v_oxvgcm,z_oxvgcm
      real,dimension(nlat,oxvgcm_nlev) ::        ! zonal means
     |  tzm_oxvgcm,uzm_oxvgcm,vzm_oxvgcm,zzm_oxvgcm
!     
! Percentage of oxvgcm data to nudge with from lbc to oxvgcm_nlev
      real :: oxvgcm_nudge(mxnlev_oxvgcm)
!
      contains
!-----------------------------------------------------------------------
      subroutine oxvgcm_bndry(modeltime)
!     
! Arg:
      integer,intent(in) :: modeltime(4)
!
! Local:
      integer :: k
      logical,save :: isread(mxtime)=.false.
      integer,save :: itime0=1
      real :: zp(nlevp1)
!
! Open the file and validate dimensions:
      if (ncid==0) call readfile(oxvgcm_ncfile)
! 
! Calculate nudging factor once per run:
      if (istep==1) then
        do k=1,nlevp1
          zp(k) = zibot+(k-1)*dlev
        enddo
        oxvgcm_nudge = 0.
        do k=1,oxvgcm_nlev
          oxvgcm_nudge(k) = 
     |      cos(pi/2.*((zp(k)-zp(1))/(zp(oxvgcm_nlev)-zp(1))))**2
        enddo 
        write(6,"('oxvgcm_bndry: oxvgcm_nlev=',i4,' oxvgcm_nudge=',
     |    /,(5f15.5))") oxvgcm_nlev,oxvgcm_nudge
      endif
!
! Read data for current time:
      if (.not.isread(itime0)) then
        call readdata(ncid,itime0)
        isread(itime0) = .true.
      endif

      end subroutine oxvgcm_bndry
!-----------------------------------------------------------------------
      subroutine readfile(ncfile)
!
! Args:
      character(len=*),intent(in) :: ncfile
!
! Local:
      integer :: istat,id,idunlim
      character(len=1024) :: msg
      integer :: nlat_rd,nlon_rd,nlev_rd
      character(len=NF90_MAX_NAME) :: dimname
!
      istat = nf90_open(ncfile,NF90_NOWRITE,ncid)
      if (istat /= NF90_NOERR) then
        write(msg,"('Error opening oxvgcm_ncfile ',a)") 
     |    trim(oxvgcm_ncfile)
        call handle_ncerr(istat,msg)
      else
        write(6,"('Opened oxvgcm_ncfile file ',a)") ncfile
      endif
!
! Get and verify latitude dimension.
! (read nlat_rd, compare with nlat, which is in params.h):
      istat = nf90_inq_dimid(ncid,'lat',id)
      istat = nf90_inquire_dimension(ncid,id,dimname,nlat_rd)
      if (nlat_rd /= nlat) then
        write(6,"(/,'>>> readfile: bad nlat_rd=',i3,
     |    ' -- should be nlat=',i3)") nlat_rd,nlat
        call shutdown('oxvgcm readfile')
      endif
!
! Get and verify longitude dimension.
! (read nlon_rd, compare with nlon, which is in params.h):
      istat = nf90_inq_dimid(ncid,'lon',id)
      istat = nf90_inquire_dimension(ncid,id,dimname,nlon_rd)
      if (nlon_rd /= nlon) then
        write(6,"(/,'>>> readfile: bad nlon_rd=',i3,
     |    ' -- should be nlon=',i3)") nlon_rd,nlon
        call shutdown('oxvgcm readfile')
      endif
!
! Get and verify levels dimension.
! (read nlev_rd, compare with mxnlev_oxvgcm, which is a parameter ',
!  in the params module)
!
      istat = nf90_inq_dimid(ncid,'lev',id)
      istat = nf90_inquire_dimension(ncid,id,dimname,nlev_rd)
      if (nlev_rd /= mxnlev_oxvgcm) then
        write(6,"(/,'>>> readfile: mxnlev_oxvgcm must be same as ',
     |    'the lev dimension on the oxvgcm file.')")
        write(6,"(  '    nlev_rd=',i4,' mxnlev_oxvgcm=',i4)")
     |    nlev_rd,mxnlev_oxvgcm
        call shutdown('oxvgcm readfile')
      endif
!
! Get number of times (value of unlimited dimension) on the file:
      istat = nf90_inq_dimid(ncid,'time',idunlim) ! id of unlimited record var
      istat = nf90_inquire_dimension(ncid,idunlim,dimname,ntime)
      if (ntime > mxtime) then
        write(6,"('>>> Reading oxvgcm data file: ntime=',i4,' > ',
     |    'mxtime=',i4)") ntime,mxtime
        call shutdown('oxvgcm ntime')
      endif

      end subroutine readfile
!-----------------------------------------------------------------------
      subroutine readdata(ncid,itime)
!
!  swb:  6/18/14:  t offset added to better reflect real atmosphere
!                  near 68-70 km  (VIRA thusfar)
! Args:
      integer,intent(in) :: ncid,itime
!
! Local:
      integer :: istat,id,i,ii
      real,dimension(nlon,nlat,mxnlev_oxvgcm) :: t,u,v,z

      istat=nf90_inq_varid(ncid,"t_lbc",id)
      istat=nf90_get_var(ncid,id,t,(/1,1,1,itime/),
     |  (/nlon,nlat,mxnlev_oxvgcm,1/))
      istat=nf90_inq_varid(ncid,"u_lbc",id)
      istat=nf90_get_var(ncid,id,u,(/1,1,1,itime/),
     |  (/nlon,nlat,mxnlev_oxvgcm,1/))
      istat=nf90_inq_varid(ncid,"v_lbc",id)
      istat=nf90_get_var(ncid,id,v,(/1,1,1,itime/),
     |  (/nlon,nlat,mxnlev_oxvgcm,1/))
      istat=nf90_inq_varid(ncid,"z_lbc",id)
      istat=nf90_get_var(ncid,id,z,(/1,1,1,itime/),
     |  (/nlon,nlat,mxnlev_oxvgcm,1/))

!     t_oxvgcm(3:nlon+2,:,:) = t             ! already in deg K
      t_oxvgcm(3:nlon+2,:,:) = t + 13.0      ! already in deg K (offset)
      u_oxvgcm(3:nlon+2,:,:) = u*100.        ! m/s to cm/s
      v_oxvgcm(3:nlon+2,:,:) = v*100.        ! m/s to cm/s
      z_oxvgcm(3:nlon+2,:,:) = z*100.        ! m to cm
!
! Periodic points:
      do i=1,2
        t_oxvgcm(i,:,:) = t_oxvgcm(nlon+i,:,:)     ! i=1,2 <- 73,74
        u_oxvgcm(i,:,:) = u_oxvgcm(nlon+i,:,:)     ! i=1,2 <- 73,74
        v_oxvgcm(i,:,:) = v_oxvgcm(nlon+i,:,:)     ! i=1,2 <- 73,74
        z_oxvgcm(i,:,:) = z_oxvgcm(nlon+i,:,:)     ! i=1,2 <- 73,74

        t_oxvgcm(nlon+2+i,:,:) = t_oxvgcm(i+2,:,:) ! i=75,76 <- 3,4
        u_oxvgcm(nlon+2+i,:,:) = u_oxvgcm(i+2,:,:) ! i=75,76 <- 3,4
        v_oxvgcm(nlon+2+i,:,:) = v_oxvgcm(i+2,:,:) ! i=75,76 <- 3,4
        z_oxvgcm(nlon+2+i,:,:) = z_oxvgcm(i+2,:,:) ! i=75,76 <- 3,4
      enddo

      if (oxvgcm_zm > 0) call calc_zm

      write(6,"('oxvgcm readdata: itime=',i4)") itime
      write(6,"('  t_oxvgcm min,max=',2e12.4)") 
     |  minval(t_oxvgcm),maxval(t_oxvgcm)
      write(6,"('  u_oxvgcm min,max=',2e12.4)") 
     |  minval(u_oxvgcm),maxval(u_oxvgcm)
      write(6,"('  v_oxvgcm min,max=',2e12.4)") 
     |  minval(v_oxvgcm),maxval(v_oxvgcm)
      write(6,"('  z_oxvgcm min,max=',2e12.4)") 
     |  minval(z_oxvgcm),maxval(z_oxvgcm)

      end subroutine readdata
!-----------------------------------------------------------------------
      subroutine calc_zm
!
! Calculate zonal means of t,u,v,z_oxvgcm, returning tzm,uzm,vzm,zzm_oxvgcm
! These arrays are in module data above.
!
! Local:
      integer :: i,j,k
      real :: rlon
      real :: tzm_diag(nlevp1,nlonp4,nlat)
      real :: uzm_diag(nlevp1,nlonp4,nlat)
      real :: vzm_diag(nlevp1,nlonp4,nlat)

      rlon = 1./real(nlon) ! note 3->nlon+2 = nlon
      tzm_oxvgcm=0. ; uzm_oxvgcm=0. ; vzm_oxvgcm=0. ; zzm_oxvgcm=0.
      do k=1,oxvgcm_nlev
!
! Output zonal means of interpolated data:
        do j=1,nlat
          do i=3,nlon+2 ! data only (not periodic points)
            tzm_oxvgcm(j,k) = tzm_oxvgcm(j,k)+t_oxvgcm(i,j,k)
            uzm_oxvgcm(j,k) = uzm_oxvgcm(j,k)+u_oxvgcm(i,j,k)
            vzm_oxvgcm(j,k) = vzm_oxvgcm(j,k)+v_oxvgcm(i,j,k)
            zzm_oxvgcm(j,k) = zzm_oxvgcm(j,k)+z_oxvgcm(i,j,k)
          enddo ! i=3,nlon+2
          tzm_oxvgcm(j,k) = tzm_oxvgcm(j,k)*rlon
          uzm_oxvgcm(j,k) = uzm_oxvgcm(j,k)*rlon
          vzm_oxvgcm(j,k) = vzm_oxvgcm(j,k)*rlon
          zzm_oxvgcm(j,k) = zzm_oxvgcm(j,k)*rlon
        enddo ! j=1,nlat
      enddo ! k=1,oxvgcm_nlev
!
! Save zm fields to secondary history:
!     do j=1,nlat
!       do i=3,nlon+2
!         tzm_diag(:,i,j) = tzm_oxvgcm(j,:) ! redundant in lon
!         uzm_diag(:,i,j) = uzm_oxvgcm(j,:) ! redundant in lon
!         vzm_diag(:,i,j) = vzm_oxvgcm(j,:) ! redundant in lon
!       enddo

!       tzm_diag(:,1:2,j) = tzm_diag(:,nlon+1:nlon+2,j)
!       tzm_diag(:,nlon+3:nlon+4,j) = tzm_diag(:,3:4,j)
!       uzm_diag(:,1:2,j) = uzm_diag(:,nlon+1:nlon+2,j)
!       uzm_diag(:,nlon+3:nlon+4,j) = uzm_diag(:,3:4,j)
!       vzm_diag(:,1:2,j) = vzm_diag(:,nlon+1:nlon+2,j)
!       vzm_diag(:,nlon+3:nlon+4,j) = vzm_diag(:,3:4,j)

!       call addfld('tzm_diag',' ',' ',tzm_diag(:,:,j),
!    |    'lev',1,nlevp1,'lon',1,nlonp4,j)
!       call addfld('uzm_diag',' ',' ',uzm_diag(:,:,j),
!    |    'lev',1,nlevp1,'lon',1,nlonp4,j)
!       call addfld('vzm_diag',' ',' ',vzm_diag(:,:,j),
!    |    'lev',1,nlevp1,'lon',1,nlonp4,j)
!     enddo
!
      end subroutine calc_zm
!-----------------------------------------------------------------------
      subroutine handle_ncerr(istat,msg)
!
! Handle a netcdf lib error:
!
      integer,intent(in) :: istat
      character(len=*),intent(in) :: msg
!
      write(6,"(/72('-'))")
      write(6,"('>>> Error from netcdf library:')")
      write(6,"(a)") trim(msg)
      write(6,"('istat=',i5)") istat
      write(6,"(a)") nf90_strerror(istat)
      write(6,"(72('-')/)")
      return
      end subroutine handle_ncerr
!-----------------------------------------------------------------------
      end module oxvgcm
