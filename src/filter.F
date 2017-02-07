!
      module filter_module
      use params_module,only: nlon,nlonp1,nlonp2,nlonp4,nlat
      implicit none
!
! Coefficients and factors for fft.
! Setfft is called once per run from init (see init_module.F)
!   Setfft calls set99, which returns trigs and ifax. 
!   Trigs and ifax are used by fftrans (see util.F) calls by filter.f 
!     and filter2.f.
!
      integer,parameter :: ntrigs = 3*nlon/2+1
      real :: trigs(ntrigs) ! e.g., if nlon==72, ntrigs==109
      integer :: ifax(13)
!
      contains
!-----------------------------------------------------------------------
      subroutine filter(f,lev0,lev1,kutj,lat)
!
! Remove longitudinal waves of prognostic variables with global fft.
! Remove wave numbers greater than kutj (see kut(nlat) in cons.F)
! This is called after mp_gatherlons, and only by tasks with mytidi==0.
! On entry, task must have global longitude data defined (mp_gatherlons).
!
! Args:
      integer,intent(in) :: kutj,lev0,lev1,lat
      real,intent(inout) :: f(nlonp4,lev0:lev1)
!
! Local:
      integer :: n1,n2,k,kk,i,ii,nlevs,nx,nw
      real :: fx(nlonp4,2*(lev1-lev0+1)),
     |        wfft((nlonp4+1)*2*(lev1-lev0+1))
!
      nlevs = lev1-lev0+1
      n1 = 2*kutj+3 ! nyquist freq (?)
      n2 = nlon+2
      if (n1 > n2) then
!       write(6,"('filter: lat=',i2,' kutj=',i2,' n1,2=',
!    |    2i3,' n1 > n2')") lat,kutj,n1,n2
        return
      endif
!
! Load fx from f for the fft:
      fx(:,:) = 0.
      do k=lev0,lev1
        do i=1,nlon
          fx(i,k) = f(i+2,k)
        enddo
      enddo
!
! Forward transform gridpoint to fourier:
! (fftrans is in util.F)
      nx = nlonp4*(2*(lev1-lev0+1))
      nw = (nlonp4+1)*2*(lev1-lev0+1)
      call fftrans(fx,nx,wfft,nw,trigs,ntrigs,ifax,1,nlonp4,nlon,nlevs,
     |  -1)
!
! Remove wave numbers greater than kutj
      do k = 1,nlevs
        do i=n1,n2
          fx(i,k) = 0.0
        enddo
      enddo
!
! Inverse transform fourier back to gridpoint:
      call fftrans(fx,nx,wfft,nw,trigs,ntrigs,ifax,1,nlonp4,nlon,nlevs,
     |  1)
!
! Redefine f from fx:
      do k=lev0,lev1
        do i=1,nlon
          f(i+2,k) = fx(i,k)
        enddo
      enddo
      end subroutine filter
!-----------------------------------------------------------------------
#include <defs.h>
      subroutine filter2(f,lev0,lev1,kutj,lat)
      use cons_module,only: dlamda
      use params_module,only: nlon,dlon
!
! Args:
      integer,intent(in) :: kutj,lev0,lev1,lat
      real,intent(inout) :: f(nlonp4,lev0:lev1)
!
! Local:
      integer,parameter :: nlond4 = nlon/4, nlond4_5 = 72/4
      integer :: nn(nlond4),nlevs,k,i,nx,nw
      real :: fx(nlonp4,2*(lev1-lev0+1)),
     |        wfft((nlonp4+1)*2*(lev1-lev0+1))
!
#if (NLAT==36 && NLON==72)
      nn=(/45,20,11,7,5,4,3,2,2,1,1,1,1,1,1,1,1,1/) ! 5.0 deg resolution
#elif (NLAT==72 && NLON==144)
      nn=(/90,90,40,40,22,22,14,14,10,10,8,8,6,6,4,4,2,2,
     |      1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1,1,1,1,1,1,1,1/) ! 2.5 deg resolution
#elif (NLAT==18 && NLON==72)              /* vtgcm single hemisphere */
      nn=(/45,20,11,7,5,4,3,2,2,1,1,1,1,1,1,1,1,1/)
#else
      UNKNOWN NLAT,NLON      ! compilation will stop here if unknown res
#endif
!
      nlevs = lev1-lev0+1
      if (kutj > nlat/2) return
      do k=lev0,lev1
        do i=1,nlon
          fx(i,k) = f(i+2,k)
        enddo
        do i=nlonp1,nlonp2
          fx(i,k)=0.
        enddo
      enddo
!
! Forward transform gridpoint to fourier:
! (fftrans is in util.F)
! subroutine fftrans(a,na,work,nw,trigs,ntrigs,ifax,inc,jump,n,lot,isign)
!
      nx = nlonp4*(2*(lev1-lev0+1))
      nw = (nlonp4+1)*2*(lev1-lev0+1)
      call fftrans(fx,nx,wfft,nw,trigs,ntrigs,ifax,1,nlonp4,nlon,nlevs,
     |  -1)
!
! Multiply by smoothing function:
      do k=lev0,lev1
        do i=1,nlon
          fx(i+2,k) = fx(i+2,k)*cos(((i+1)/2)*dlamda/2.)**(2*nn(kutj))
        enddo ! i=1,nlon
      enddo ! k=lev0,lev1
!
! Inverse transform back to gridpoint:
      call fftrans(fx,nx,wfft,nw,trigs,ntrigs,ifax,1,nlonp4,nlon,nlevs,
     |  +1)
!
! Save smoothed field:
      do k=lev0,lev1
        do i=1,nlon
          f(i+2,k) = fx(i,k)
        enddo ! i=1,nlon
      enddo ! k=lev0,lev1
      end subroutine filter2
!-----------------------------------------------------------------------
      subroutine filter_sp(fout,lev0,lev1,lon0,lon1,lat0,lat1,kut)
!
! Filter a composition species:
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
!     code = 123 ; state = 'filter_sp' ; activity='Filtering'
      call vtbegin(123,ier)
#endif
!     write(6,"('Enter filter_sp')")
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
!     code = 123 ; state = 'filter_sp' ; activity='Filtering'
      call vtend(123,ier)
#endif
      end subroutine filter_sp
!-----------------------------------------------------------------------
      end module filter_module
