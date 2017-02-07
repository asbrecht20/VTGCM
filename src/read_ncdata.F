
      module read_ncdata_module
      use nchist_module,only: nc_open,nc_close,handle_ncerr
      use params_module,only: nlat,nlon,nlonp2,nlonp4,nlevp1
      implicit none
#include <netcdf.inc>
!
      real,dimension(nlonp4,nlat,nlevp1) :: q,ui,vi,rj,quv
      real,dimension(nlevp1) :: tm1,tm2,xnod,xnco2d,xnon,xnco2n,
     |  cool1,cool2
!
! Fields derived from above vars, for use in sub cooling (dt.F):
      real,dimension(nlevp1) :: 
     |  tm1a,tm2a, ! tm1,tm2 shifted 1/2 dz             ! T5,T6
     |  aday,      ! formerly AD(ZKMXP)
     |  anight     ! formerly AN(ZKMXP)
!
! Change cooling constant, as per Bougher, 9/29/06:
!     real,parameter :: coolk = 1.1550e-13  ! old value
      real,parameter :: coolk = 1.7321e-13  ! new value
!
      contains
!-----------------------------------------------------------------------
      subroutine read_ncdata(ncfile)
      implicit none
!
! Read data volume.
!
! Args:
      character(len=*),intent(in) :: ncfile
!
! Local:
      integer,parameter :: mxdims=20
      integer :: ncid,istat,ndims,nvars,ngatts,idunlim,i,natts,itype,
     |  iddims(mxdims),idlev,nlevrd
      character(len=120) :: diskfile
      character(len=80) :: varname
      real :: rdlon(nlon),rdlat(nlat),rdlev(nlevp1)
      real,dimension(nlon,nlat,nlevp1) :: rdvar ! local var for read
!
      write(6,"('read_ncdata: ncfile=',a)") trim(ncfile)
      diskfile = ' '
      call getfile(ncfile,diskfile)

      call nc_open(ncid,diskfile,'OLD','READ')
      write(6,"('read_ncdat: opened file ',a,' ncid=',i3)")
     |  trim(diskfile),ncid

      istat = nf_inq(ncid,ndims,nvars,ngatts,idunlim)

      write(6,"('read_ncdata: ndims=',i3,' nvars=',i3,' ngatts=',i3,
     |  ' idunlim=',i3)") ndims,nvars,ngatts,idunlim

      istat = nf_inq_dimid(ncid,'lev',idlev)
      istat = nf_inq_dimlen(ncid,idlev,nlevrd)
      if (nlevrd /= nlevp1) then
        write(6,"('>>> read_ncdata: nlevrd=',i4,' but nlevp1=',i4)")
     |    nlevrd,nlevp1
        call shutdown('read_ncdata') 
      endif

      do i=1,nvars
        istat = nf_inq_var(ncid,i,varname,itype,ndims,iddims,natts)
        if (istat /= NF_NOERR) then
          write(6,"('>>> nc_rdhist: error ',
     |      'inquiring about var: i=',i4,' nvars=',i4)") i,nvars
          call handle_ncerr(istat,'inquiring about a var')
        endif
!       write(6,"(/,'read_ncdata: i=',i3,' field=',a,' itype=',i3,
!    |    ' ndims=',i3,' natts=',i2,' iddims=',/,20i3)") 
!    |    i,trim(varname),itype,ndims,natts,iddims
        select case (trim(varname))
          case('lon')
!           istat = nf_get_var
          case('lat')
          case('lev')
!           call rdf1d(ncid,i,rdlev,'lev')
!
! 3d reads (nlon,nlat,nlevp1), and transfer to module data (nlonp4,nlat,nlevp1):
          case('Q') 
            call rdf3d(ncid,i,rdvar,'Q')
            q(3:nlonp2,:,:) = rdvar(:,:,:)
          case('UI') 
            call rdf3d(ncid,i,rdvar,'UI')
            ui(3:nlonp2,:,:) = rdvar(:,:,:)
          case('VI') 
            call rdf3d(ncid,i,rdvar,'VI')
            vi(3:nlonp2,:,:) = rdvar(:,:,:)
          case('RJ') 
            call rdf3d(ncid,i,rdvar,'RJ')
            rj(3:nlonp2,:,:) = rdvar(:,:,:)
          case('QUV') 
            call rdf3d(ncid,i,rdvar,'QUV')
            quv(3:nlonp2,:,:) = rdvar(:,:,:)
!
! 1-d column reads:
          case('TM1') 
            call rdf1d(ncid,i,tm1,'TM1')
          case('TM2') 
            call rdf1d(ncid,i,tm2,'TM2')
          case('XNOD') 
            call rdf1d(ncid,i,xnod,'XNOD')
          case('XNCO2D') 
            call rdf1d(ncid,i,xnco2d,'XNCO2D')
          case('XNON') 
            call rdf1d(ncid,i,xnon,'XNON')
          case('XNCO2N') 
            call rdf1d(ncid,i,xnco2n,'XNCO2N')
          case('COOL1') 
            call rdf1d(ncid,i,cool1,'COOL1')
          case('COOL2') 
            call rdf1d(ncid,i,cool2,'COOL2')
          case default
            write(6,"('read_ncdata: unused var: ',a)") trim(varname)
        end select
      enddo
      call nc_close(ncid)
!
! Periodic points for 3d fields: 
      q(1:2,:,:) = q(nlon+1:nlonp2,:,:)
      q(nlon+3:nlonp4,:,:) = q(3:4,:,:)
      ui(1:2,:,:) = ui(nlon+1:nlonp2,:,:)
      ui(nlon+3:nlonp4,:,:) = ui(3:4,:,:)
      vi(1:2,:,:) = vi(nlon+1:nlonp2,:,:)
      vi(nlon+3:nlonp4,:,:) = vi(3:4,:,:)
      rj(1:2,:,:) = rj(nlon+1:nlonp2,:,:)
      rj(nlon+3:nlonp4,:,:) = rj(3:4,:,:)
      quv(1:2,:,:) = quv(nlon+1:nlonp2,:,:)
      quv(nlon+3:nlonp4,:,:) = quv(3:4,:,:)
!
! Calculate aday, anight:
      call lincoeff
!
      end subroutine read_ncdata
!-----------------------------------------------------------------------
      subroutine rdf3d(ncid,id,f,name)
      implicit none
!
! Args:
      integer,intent(in) :: ncid,id
      real,dimension(nlon,nlat,nlevp1),intent(out) :: f
      character(len=*),intent(in) :: name 
!
! Local:
      integer :: istat,idlocal
!
      idlocal = id
      istat = nf_get_var_double(ncid,idlocal,f)
      if (istat /= NF_NOERR) then
        write(6,"('>>> rdf3d: error reading 3d var ',a,': id=',i4)") 
     |    trim(name),id
        call handle_ncerr(istat,'reading 3d var')
      else
        write(6,"('rdf3d: read field ',a,' id=',i3,' min,max=',2e12.4)") 
     |    trim(name),id,minval(f),maxval(f)
      endif
      end subroutine rdf3d
!-----------------------------------------------------------------------
      subroutine rdf1d(ncid,id,f,name)
      implicit none
!
! Args:
      integer,intent(in) :: ncid,id
      real,dimension(nlevp1) :: f
      character(len=*),intent(in) :: name 
!
! Local:
      integer :: istat
!
      istat = nf_get_var_double(ncid,id,f)
      if (istat /= NF_NOERR) then
        write(6,"('>>> rdf1d: error reading 1d var ',a,': id=',i4)") 
     |    trim(name),id
        call handle_ncerr(istat,'reading 1d var')
      else
!       write(6,"('rdf1d: read field ',a,' = ',/,(6e12.4))")
!    |    trim(name),f
      endif

      end subroutine rdf1d
!-----------------------------------------------------------------------
      subroutine lincoeff
      use params_module,only: dz,zibot
      implicit none
!
! Local:
      integer :: k
      real :: roco2,roco2m,xxx,efd,efn,gee,efnm,facm,fden,stm,drate2
      real :: szp(nlevp1)
!
! Shift tm1,2 by 1/2 level:
      do k=1,nlevp1-1
        tm1a(k) = .5*(tm1(k)+tm1(k+1))
        tm2a(k) = .5*(tm2(k)+tm2(k+1))
      enddo
      tm1a(nlevp1) = tm1a(nlevp1-1)
      tm2a(nlevp1) = tm2a(nlevp1-1)
      do k=1,nlevp1
        szp(k) = zibot + .5*dz + (k-1)*dz
      enddo
!
! 9/20/06 btf: The next 2 k-loops are taken from old vtgcm code
!              (/home/tgcm/vtgcm/src.1996/start.f) 
! 10/31/06 swb:  Checked new code with vtgcm2 code. Looks fine.
!
C     ***************************************************************
C     DAYSIDE LIN COEFF CALC
C     ***************************************************************
      do k=1,nlevp1                                  ! DO 261 in start.f
        roco2m = coolk*sqrt(tm1a(k))/(5.2e-15*(1.+2.e-3*
     |    (tm1a(k)-220.))*(tm1a(k)/273.3))
        efd = exp(szp(k)+1.5)*((xnco2d(k)+xnod(k))/(xnco2d(k)+xnod(k)*
     |    roco2m))
        efnm = efd
        facm = (1.+efd)/(1.+efnm)
        fden = ((1.+2.e-3*(tm1a(k)-220.))*(tm1a(k)/273.3))
        stm = sqrt(tm1a(k))
        drate2 = (coolk/5.2e-15)*(1./(2.*stm*fden)-(stm/273.3)*
     |    (1.+2.e-3*(2.*tm1a(k)-220.))/fden**2.)
        aday(k) = facm*cool1(k)*((960./tm1a(k)**2)-     ! AD(K))
     |    (efnm/(1.+efnm))*(40.6/(3.*tm1a(k)**4./3.)+
     |    xnod(k)/(xnco2d(k)+xnod(k)*roco2m)*drate2))
      enddo ! k=1,nlevp1
      write(6,"('lincoeff: aday=',/,(6e12.4))") aday
C     ***************************************************************
C     NIGHTSIDE LIN C0EFF
C     ***************************************************************
      do k=1,nlevp1                                  ! DO 263 in start.f
        roco2m = coolk*sqrt(tm2a(k))/(5.2e-15*(1.+2.e-3*
     |    (tm2a(k)-220.))*(tm2a(k)/273.3))
        efd = exp(szp(k)+1.5)*((xnco2n(k)+xnon(k))/(xnco2n(k)+xnon(k)*
     |    roco2m))
        efnm = efd
        facm = (1.+efd)/(1.+efnm)
        fden = ((1.+2.e-3*(tm2a(k)-220.))*(tm2a(k)/273.3))
        stm = sqrt(tm2a(k))
        drate2 = (coolk/5.2e-15)*(1./(2.*stm*fden)-(stm/273.3)*
     |    (1.+2.e-3*(2.*tm2a(k)-220.))/fden**2.)
        anight(k) = facm*cool2(k)*((960./tm2a(k)**2)-   ! AN(K))
     |    (efnm/(1.+efnm))*(40.6/(3.*tm2a(k)**4./3.)+
     |    xnon(k)/(xnco2n(k)+xnon(k)*roco2m)*drate2))
      enddo ! k=1,nlevp1
      write(6,"('lincoeff: anight=',/,(6e12.4))") anight
      end subroutine lincoeff
!-----------------------------------------------------------------------
      end module read_ncdata_module
