!
      subroutine allocdata
!
! Allocate msc data. Main fields are not allocated here, see
! init_fields in fields.F. This is not in field.F in order to 
! avoid circular module dependencies.
!
! 1/27/06 btf:  call alloc_mgw (mgw.F) for zp-10 lbc 
! 2/15/07 swb:  call alloc_o2 
! 03/12/08 swb: call alloc_n4s 
! 05/12/08 swb: call alloc_no 
! btf 6/13/13: call alloc_so, alloc_so2
!
      use params_module,only: nlevp1,nlonp4,nlat
      use fields_module,only: f4d,f3d,nf4d,nf3d,field_initval,foutput,
     |  fzg
      use mpi_module,only: lon0,lon1,lat0,lat1,mytid
      use chemrates_module,only: alloc_tdep
      use qrj_module,only: alloc_q
      use o2_module,only: alloc_o2
      use n4s_module,only: alloc_n4s, alloc_no
      use so_module,only: alloc_so
      use so2_module,only: alloc_so2
!
!     use efield_module,only: alloc_e
!     use no_module,only: alloc_no
!     use magpres_g_module,only: alloc_pg
!     use mgw,only: alloc_mgw
!     use tuv_bndry,only: alloc_tuvbnd
!     use msis_module,only: alloc_msis
      implicit none
!
! Local:
      integer :: i,istat,nlons,nlats,len 
!
! Allocate temperature-dependent chemical reaction rates:
      call alloc_tdep
!
! Allocate ionization rates (qrj module):
      call alloc_q(lon0,lon1,lat0,lat1)
!
! Allocate boundary conditions and production and loss terms
! for minor composition species:
! ** Allocate o2_ubc, o2_lbc, o2_prod, o2_loss (o2_module):
      call alloc_o2(lon0,lon1,lat0,lat1)
! ** Allocate no_ubc, no_lbc, no_prod, no_loss (n4s_module):
      call alloc_no(lon0,lon1,lat0,lat1)
! ** Allocate n4s_ubc, n4s_lbc, n4s_prod, n4s_loss (n4s_module):
      call alloc_n4s(lon0,lon1,lat0,lat1)
! ** Allocate so_ubc, so_lbc, so_prod, so_loss (so_module):
      call alloc_so(lon0,lon1,lat0,lat1)
! ** Allocate so2_ubc, so2_lbc, so2_prod, so2_loss (so2_module):
      call alloc_so2(lon0,lon1,lat0,lat1)
!
! Allocate J(mag.pressure,gravity (magpres_g module):
!     call alloc_pg(lon0,lon1,lat0,lat1)
!
! Allocate electric field (efield module):
!     call alloc_e(lon0,lon1,lat0,lat1)
!
! Allocate arrays for gravity wave parameterization:
!     call alloc_mgw(nlevp1,lon0,lon1,lat0,lat1)
!
! Allocate t,u,v lower boundary:
!     call alloc_tuvbnd(lon0,lon1,lat0,lat1)
!
! Allocate msis_tn lower boundary:
!     call alloc_msis(lon0,lon1,lat0,lat1)
!
! We are not doing parallel i/o to netcdf output histories, so when
! it comes time to write a history, the task subdomains are gathered
! to the root task (see call mp_gather2root in advance), and only
! the root task writes to the file. Here we allocate the full 3d grid
! for all history fields on the root task. This is where the data
! will be gathered and written by the root when a history is written. 
!
! The pointer is declared in fields.F:
!     real,pointer :: foutput(:,:,:,:) ! (nlevp1,nlonp4,nlat,nf4d)
! 
#ifdef MPI
      if (mytid==0) then
        allocate(foutput(nlevp1,nlonp4,nlat,nf4d),stat=istat)
        if (istat /= 0) then
          write(6,"('>>> allocdata: error allocating 3d data for ',
     |      'all fields on root task: nlevp1=',i3,' nlonp4=',i3,
     |      ' nlat=',i3,' nf4d=',i3,' istat=',i3)")
     |      nlevp1,nlonp4,nlat,nf4d,istat
        else
          write(6,"('Allocated 3d data for all fields on root task:')")
          len = nlevp1*nlonp4*nlat*nf4d
          write(6,"('  nlevp1*nlonp4*nlat*nf4d = ',i3,'*',i3,'*',
     |      i3,'*',i3,' = ',i8,' words *8 = ',i10,' bytes.')")
     |      nlevp1,nlonp4,nlat,nf4d,len,len*8
        endif
        allocate(fzg(nlevp1,nlonp4,nlat),stat=istat)
        if (istat /= 0) 
     |    write(6,"('>>> allocdata: error allocating fzg')")
      endif
#else
      allocate(foutput(nlevp1,nlonp4,nlat,nf4d),stat=istat)
      if (istat /= 0) then
        write(6,"('>>> allocdata: error allocating 3d data for ',
     |    'all fields on root task: nlevp1=',i3,' nlonp4=',i3,
     |    ' nlat=',i3,' nf4d=',i3,' istat=',i3)")
     |    nlevp1,nlonp4,nlat,nf4d,istat
      else
        write(6,"('Allocated 3d data for all fields on root task:')")
        len = nlevp1*nlonp4*nlat*nf4d
        write(6,"('  nlevp1*nlonp4*nlat*nf4d = ',i3,'*',i3,'*',
     |    i3,'*',i3,' = ',i8,' words *8 = ',i10,' bytes.')")
     |    nlevp1,nlonp4,nlat,nf4d,len,len*8
      endif
      allocate(fzg(nlevp1,nlonp4,nlat),stat=istat)
      if (istat /= 0) 
     |  write(6,"('>>> allocdata: error allocating fzg')")
#endif
      end subroutine allocdata
