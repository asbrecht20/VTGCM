pro mknc_lbc
;
; Read text files w/ lbc for vtgcm lower boundary of t,u,v,z,
; and write a netcdf file for import to vtgcm
;
;-rw-r--r-- 1 foster tgcm 272320 Nov 20 08:56 GeoHtout.dat
;-rw-r--r-- 1 foster tgcm 246934 Nov 20 08:56 Meriout.dat
;-rw-r--r-- 1 foster tgcm 220480 Nov 20 08:56 Tempout.dat
;-rw-r--r-- 1 foster tgcm 189170 Nov 20 08:56 Zonalout.dat
;
  nlon = 72
  nlat = 36
  nlev = 10 ; -16 to -11.5 by 0.5

  tfile = 'Tempout.dat'
  ufile = 'Zonalout.dat'
  vfile = 'Meriout.dat'
  zfile = 'GeoHtout.dat'

  t_lbc = fltarr(nlon,nlat,nlev)
  u_lbc = fltarr(nlon,nlat,nlev)
  v_lbc = fltarr(nlon,nlat,nlev)
  z_lbc = fltarr(nlon,nlat,nlev)

  fglb = fltarr(nlon,nlat)
;
; Levels coordinate (these are not easily read from the 
; ascii files, so hardwire them here):
;
  lev = fltarr(nlev)
  dlev = 0.5
  lev0 = -16.
  for k=0,nlev-1 do lev[k] = lev0+k*dlev
  print,'nlev=',nlev,' lev=' & print,lev
;
; Read temperature:
;
  openr,lu,tfile,/get_lun
  print,'Opened tfile ',tfile,' lu=',lu
  line = ''
  for k=0,nlev-1 do begin
    readf,lu,line
    readf,lu,fglb
    t_lbc[*,*,k] = fglb
    print,'k=',k,' t_lbc min,max=',min(t_lbc[*,*,k]),max(t_lbc[*,*,k])
  endfor
  free_lun,lu
;
; Read zonal wind:
;
  openr,lu,ufile,/get_lun
  print,'Opened ufile ',ufile,' lu=',lu
  for k=0,nlev-1 do begin
    readf,lu,line
    readf,lu,fglb
    u_lbc[*,*,k] = fglb
    print,'k=',k,' u_lbc min,max=',min(u_lbc[*,*,k]),max(u_lbc[*,*,k])
  endfor
  free_lun,lu
;
; Read meridional wind:
;
  openr,lu,vfile,/get_lun
  print,'Opened vfile ',vfile,' lu=',lu
  for k=0,nlev-1 do begin
    readf,lu,line
    readf,lu,fglb
    v_lbc[*,*,k] = fglb
    print,'k=',k,' v_lbc min,max=',min(v_lbc[*,*,k]),max(v_lbc[*,*,k])
  endfor
  free_lun,lu
;
; Read geopotential:
;
  openr,lu,zfile,/get_lun
  print,'Opened ufile ',zfile,' lu=',lu
  for k=0,nlev-1 do begin
    readf,lu,line
    readf,lu,fglb
    z_lbc[*,*,k] = fglb
    print,'k=',k,' z_lbc min,max=',min(z_lbc[*,*,k]),max(z_lbc[*,*,k])
  endfor
  free_lun,lu
;
; Create new netcdf dataset:
;
  ncfile_name = 'oxvgcm_lbc.nc'            ; Data is from Oxford Venus GCM
  ncid = ncdf_create(ncfile_name,/clobber)
  print,'Created netcdf file ',ncfile_name,' ncid=',ncid
;
; Define dimensions:
;
  id_unlim = ncdf_dimdef(ncid,'time',/unlimited)
  id_nlon = ncdf_dimdef(ncid,'lon',nlon)
  id_nlat = ncdf_dimdef(ncid,'lat',nlat)
  id_nlev = ncdf_dimdef(ncid,'lev',nlev)
;
; Define coordinate and data variables:
;
  idv_time = ncdf_vardef(ncid,'time',id_unlim)
  idv_lon  = ncdf_vardef(ncid,'lon' ,id_nlon)
  idv_lat  = ncdf_vardef(ncid,'lat' ,id_nlat)
  idv_lev  = ncdf_vardef(ncid,'lev' ,id_nlev)

  idv_t    = ncdf_vardef(ncid,'t_lbc',[id_nlon,id_nlat,id_nlev,id_unlim])
  idv_u    = ncdf_vardef(ncid,'u_lbc',[id_nlon,id_nlat,id_nlev,id_unlim])
  idv_v    = ncdf_vardef(ncid,'v_lbc',[id_nlon,id_nlat,id_nlev,id_unlim])
  idv_z    = ncdf_vardef(ncid,'z_lbc',[id_nlon,id_nlat,id_nlev,id_unlim])
;
; Take out of define mode:
;
  ncdf_control,ncid,/endef
;
; Write coordinate variables:
;
; Time coordinate:
  time = 1.                      ; arbitrary for now
  ncdf_varput,ncid,idv_time,time
;
; Longitude coordinates:
;
  lon = fltarr(nlon)
  dlon = 5.0
  for i=0,nlon-1 do begin
    lon[i] = -180.+i*dlon
  endfor
  print,'lon=' & print,lon
  ncdf_varput,ncid,idv_lon,lon
;
; Latitude coordinates:
;
  lat = fltarr(nlat)
  dlat = 5.0
  for j=0,nlat-1 do begin
    lat[j] = -87.5+j*dlat
  endfor
  print,'lat=' & print,lat
  ncdf_varput,ncid,idv_lat,lat
;
; Level coordinates (lev is set above):
;
  ncdf_varput,ncid,idv_lev,lev
;
; Write data variables:
;
  ncdf_varput,ncid,idv_t,t_lbc
  ncdf_varput,ncid,idv_u,u_lbc
  ncdf_varput,ncid,idv_v,v_lbc
  ncdf_varput,ncid,idv_z,z_lbc

  ncdf_close,ncid
end
