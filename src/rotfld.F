      module rotate_field
!
! Rotate a field on the longitude grid to follow planet rotation in
! elapsed time segment.
!
      implicit none
      integer,parameter :: 
!    |  nsecperday = 20997000 ! secs per day planet rotation (assume 24 hrs)
     |  nsecperday = 3600*24 ! secs per day planet rotation (assume 24 hrs)
      real :: deg_per_sec
      integer,private ::
     |  isecprot=0   ! secs at previous rotation
      contains
!-----------------------------------------------------------------------
      subroutine rotf3d(f,fname,lev0,lev1,lon0,lon1,lat0,lat1,step,
     |  istep,update)
!
! Rotate field f, which is on subdomains, by gathering longitudes to
! leftmost mpi tasks (left longitude column, mytidi==0), where they are 
! rotated globally by sub rotfld, then scattered back out to the other 
! task subdomains.
!
      use params_module,only: nlon,nlonp4
#ifdef MPI
      use mpi_module,only: mp_gatherlons_f3d,mp_scatterlons_f3d,mytidi
      implicit none
#else
      implicit none
      integer :: mytidi=0
#endif
!
! Args:  
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1,step,istep
      real,intent(inout) :: f(lev0:lev1,lon0:lon1,lat0:lat1)
      logical,intent(in) :: update
      character(len=*),intent(in) :: fname ! diag only
!
! Local:
      integer :: i,j,isec_curr,ier
      real :: fkij(lev0:lev1,nlonp4,lat0:lat1)
      real ::  
     |  deg_per_lon,  ! degrees rotation per longitude grid point
     |  deg_elapsed   ! degrees rotation since last rotation
!
      isec_curr = step*istep               ! seconds from beginning of run
      deg_per_lon = 360./float(nlon)       ! deg rotation per grid point
      deg_per_sec = 360./float(nsecperday) ! deg rotation per sec
!
! Elapsed rotation since last rotation update
      deg_elapsed  = deg_per_sec*(isec_curr-isecprot) 
!
      if (deg_elapsed >= deg_per_lon) then 
!
! Define lons from current task:
        fkij = 0.
        do j=lat0,lat1
          do i=lon0,lon1
            fkij(:,i,j) = f(:,i,j)
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
! Only leftmost tasks at each j-row of tasks does the rotation, because
! they have the full global longitude dimension:
!
        if (mytidi==0) then

!         if (trim(fname)=='TN')
!    |    write(6,"('call rotfld: deg_per_lon=',f8.2,' deg_elapsed=',
!    |      f8.2,' isecprot=',i8,' isec_curr=',i8,' elapsed secs=',i8)")
!    |      deg_per_lon,deg_elapsed,isecprot,isec_curr,
!    |      isec_curr-isecprot

          call rotfld(isecprot,isec_curr,fkij(lev0,1,lat0),
     |      fname,lev1-lev0+1,nlonp4,lat1-lat0+1,nlonp4,nlon,
     |      istep,ier)
          if (ier > 0) then
            write(6,"(/,'>>> rotf3d: error from rotfld')")
            call shutdown('rotf3d')
          endif
!
! Periodic points:
!         fkij(:,1:2,:) = fkij(:,nlon+1:nlon+2,:) ! 1:2 <- 73:74
!         fkij(:,nlon+3:nlon+4,:) = fkij(:,3:4,:) ! 75:76 <- 3:4
        endif ! mytidi==0
#ifdef MPI
!
! Now leftmost task at each j-row must redistribute rotated data
! back to other tasks in the j-row (mytidi>0,mytidj) (includes latitude):
!
        call mp_scatterlons_f3d(fkij,lev0,lev1,lon0,lon1,lat0,lat1,1)
#endif
!
! Return scattered and rotated array to f:
        do j=lat0,lat1
          do i=lon0,lon1
            f(:,i,j) = fkij(:,i,j)
          enddo
        enddo ! j=lat0,lat1
!
! Update time of previous rotation (module data):
        if (update) then
          isecprot = isec_curr
        endif
!
! Do not rotate because elapsed rotation is less than degrees per gridpoint:
      else ! deg_elapsed < deg_per_lon
!       write(6,"('rotf3d: Did *not* rotate: deg_elapsed=',f8.2,
!    |    ' deg_per_lon=',f8.2)") deg_elapsed,deg_per_lon 
      endif
      end subroutine rotf3d
!-----------------------------------------------------------------------
      subroutine rotfld(isec_prev,isec_curr,f,fname,id1,id2,id3,nlondim,
     |  nlon,istep,ier)
      implicit none
!
! Rotate field f(id1,id2,id3) in longitude, from previous to current
! time. nlondim is the longitude dimension (must be == id1, id2, or id3)
! Assume longitude dimension is global west to east, and rotation 
! is performed west to east, i.e., following planet rotation.
!
! WARNING: This routine assumes nlondim is the *global* longitude dimension,
!          (including periodic points, e.g., nlondim=76 for 5 deg grid).
!          It will *not* work on longitude subdomains. Field f must have 
!          data on the global longitude grid. See mpi.F for gather/scatter
!          utilities, e.g., mp_gatherlons_f3d and mp_scatterlons_f3d.
!
! Args:
      integer, intent(in) :: 
     |  isec_prev,    ! seconds at previous rotation
     |  isec_curr,    ! current seconds (elapsed = curr-prev)
     |  id1,id2,id3,  ! dimensions of f
     |  nlondim,      ! lon dimension of f (must be == id1, id2, or id3)
     |  nlon,         ! lon dimension, excluding periodic points (will be nlondim-4)
     |  istep         ! diag only
      real,intent(inout) :: f(id1,id2,id3)
      integer,intent(out) :: ier ! error has occurred if ier > 0
      character(len=*),intent(in) :: fname
!
! Local:
      integer :: i,ii,iii,nlonrot,newlon(nlondim)
      real :: deg_rot,degperlon
      real :: ftmp(id1,id2,id3)
!
! nlondim must be one of the 3 dimensions of f:
      ier = 0
      if (nlondim /= id1 .and. nlondim /= id2 .and. nlondim /= id3) then
        write(6,"('>>> rotfld: nlondim must == id1,id2, or id3: ',
     |    'nlondim=',i4,' id1,2,3=',3i5)") nlondim,id1,id2,id3
        ier = 1
        return
      endif
!
! Degrees of planetary rotation per sec:
      deg_per_sec = 360./float(nsecperday)
!
! Degrees of rotation in elapsed time since last rotation:
      deg_rot = float(isec_curr-isec_prev)*deg_per_sec
!
! nlonrot = number of longitude grid points to rotate the field:
      degperlon = 360./float(nlon) ! deg rotation per lon grid point
      nlonrot = nint(deg_rot/degperlon)
      if (nlonrot <= 0) then
        write(6,"('>>> Warning rotfld: nlonrot=',i3,
     |    ' -- no rotation was done.')") nlonrot
        ier = 1
        return
      elseif (nlonrot > nlon) then
        write(6,"('>>> Warning rotfld: nlonrot must be <= nlon:',
     |    ' nlonrot=',i4,' nlon=',i4)") nlonrot,nlon
        ier = 1
        return
      else
        if (trim(fname)=='TN')
     |  write(6,"('rotfld: field =',a,' istep=',i4,' nlonrot=',i3,
     |    ' deg_rot=',f8.4,' degperlon=',f8.4)") fname(1:8),
     |    istep,nlonrot,deg_rot,degperlon
      endif
!
! Rotate the field nlonrot longitude grid points:
      do i=1,nlondim
        newlon(i) = mod(i+nlonrot-1,nlondim)+1 ! destination lon index
      enddo
!     write(6,"('rotfld: newlon=',/,(15i4))") newlon
      if (nlondim==id1) then
        do i=1,nlondim
          ftmp(i,:,:) = f(newlon(i),:,:)
        enddo
      elseif (nlondim==id2) then ! this is typical (f4d arrays)
        do i=1,nlondim
          ftmp(:,i,:) = f(:,newlon(i),:)
        enddo
      elseif (nlondim==id3) then
        do i=1,nlondim
          ftmp(:,:,i) = f(:,:,newlon(i))
        enddo
      endif
      f = ftmp
      end subroutine rotfld
!-----------------------------------------------------------------------
      end module rotate_field
