      subroutine advance
!
! Advance the model nstep time steps.
! -- Latest revision:  swb:  05/02/08:  n2d,o2p,op,ne
!                                       to conditional for itc (OK)
! -- Latest revision:  swb:  05/02/08:  o2, o2_nm, n4s, n4s_nm
!                                       to conditional for itc (OK)
! -- Latest revision:  swb:  05/12/08:  no, no_nm, 
!                                       to conditional for itc
! -- Latest revision:  swb:  10/21/13:  so, so_nm, so2, so2_nm 
!                                       to conditional for itc
!
      use fields_module
      use hist_module,only: nstep,modeltime,nsource
      use init_module,only: istep,uthr,iter,secs,iday,iyear
      use timing_module,only: timer,timing
      use cons_module,only: dt,secperhr,secperday
      use input_module,only: start,calendar_advance,step
      use params_module,only: nlevp1,nlat,nlonp4
      use rotate_field,only: rotf3d
      use output,only: outhist
      use addfld_module,only: addfld
      use hdif_module,only: hdif1,hdif2
      use lbc,only: tuvz_lbc
#ifdef MPI
      use mpi_module,only: lat0,lat1,lon0,lon1,mytid,ntask,
     |  mp_gather2root,mp_bndlats,mp_bndlons,mp_polelats,
     |  mp_periodic_f4d
      implicit none
#include <mpif.h>
#else
      implicit none
      integer :: mytid=0,ntask=1
      integer :: lat0=1,lat1=nlat,lon0=1,lon1=nlonp4
#endif
!
! Local:
      integer :: i,j,ier,itmp,ix
      real :: fmin,fmax
      logical :: debug=.false.
      integer(kind=8) :: nsecs ! current model time in seconds,
                               ! including model day
      real :: secs_per_step,dday,
     |  time0,          time1,
     |  time0_step,     time1_step,
     |  time0_init,     time1_init,
     |  time0_phist,    time1_phist,
     |  time0_shist,    time1_shist,
     |  time0_prep,     time1_prep,
     |  time0_barrier,  time1_barrier
      character(len=3) :: timing_type
      logical :: wrprim,wrsech,save_prim,save_sech,newseries_sech,
     |  newseries_prim,iseries_sech,iseries_prim,time2write
!
! External:
      logical,external :: wrhist,time2print
!
      write(6,"(/,'Enter advance: iter=',i6,' nstep=',i6)") iter,nstep
      call timer(time0_init,time1_init,'INIT',0,0)
!
! Init timestep index, modeltime, and nsecs:
      istep = 0
      modeltime(1:3) = start(:,1)
      modeltime(4) = 0
      call modeltime_to_nsecs(modeltime,nsecs) ! sets nsecs, includes day
      timing_type = 'sys'
      if (timing%rtc) timing_type = 'rtc'
      if (debug) write(6,"('advance: nsecs=',i12)") nsecs
!
! Echo source history to primary history output if an initial run,
!   and write initial secondary history if necessary (istep==0).
! Source history will not be written if a continuation run (i.e., if
!   SOURCE was *not* provided by user input.
! (if mpi, only master writes histories)
!
! Source was read into data(:,:,:,itp), and output always writes from
! data(:,:,:,itc), so swap indices to echo itp source data:
        itmp = itp
        itp = itc
        itc = itmp
#ifdef MPI
!
! If an initial run, echo source history to output before first time step.
!   Because only subdomains are read by rdhist, and only the root task
!   outputs histories, the subdomain data must be gathered to the root.
!
        if (debug) write(6,"('advance: before mp_gather2root for ',
     |    'source')")
!
! 9/21/06 btf: code seg faults in mp_gather2root_prim on lightning.
!
        call mp_gather2root(itc,'prim')
        if (debug) write(6,"('advance after mp_gather2root for ',
     |    'source')")
        if (mytid==0) call outhist(istep,modeltime)
        if (debug) write(6,"('advance after outhist for source')")
#else
        do i=1,nf4d
          foutput(:,lon0:lon1,lat0:lat1,i) =
     |      f4d(i)%data(:,lon0:lon1,lat0:lat1,itc)
        enddo
        call outhist(istep,modeltime)
#endif
!
! nstep may be zero if user only wants to copy source history:
      if (nstep==0) then
        write(6,"(/,'ADVANCE: model is not advanced in time because ',
     |    ' start==stop.')")
        call timer(time0_init,time1_init,'INIT',1,0) ! end init timing
        return
      endif
!
! Reswap indices, so model reads itp data and updates itc data:
      itmp = itp
      itp = itc
      itc = itmp
!
      call timer(time0_init,time1_init,'INIT',1,0) ! end init timing
      if (timing%level >= 2)
     |  write(6,"('Time in INIT =          ',
     |    f12.2,' Advance: step ',i5)") time1_init,istep
!
! Main time loop:
  100 continue
#ifdef VT
!     code = 127 ; state = 'timestep' ; activity='ModelCode'
      call vtbegin(127,ier)
#endif
      call timer(time0_step,time1_step,'STEP',0,0) ! start step timing
      call timer(time0_prep,time1_prep,'PREP',0,0) ! start prep timing
      iter=iter+1
      istep = istep+1
!
! Increment model time by one time step (day,hr,min,secs):
!
      nsecs = nsecs+step
      call nsecs_to_modeltime(nsecs,modeltime) ! increments modeltime(4)
!     secs = mod(iter*int(dt),86400)
!     uthr = secs/3600.
      secs = mod(int(iter,8)*int(dt,8),int(secperday,8))
      uthr = secs/secperhr

!     write(6,"('step=',i4,' iter=',i5,' int(dt)=',i3,' modeltime=',
!    |  3i4,' secs=',f10.2,' uthr=',f8.2)") istep,iter,int(dt),
!    |  modeltime(1:3),secs,uthr

!
! Report to stdout:
      if (time2print(nstep,istep)) then
        if (istep > 1) then
          write(6,"('Step ',i6,' of ',i6,' mtime=',3i3,
     |      ' secs/step (',a,') =',f6.2)") istep,nstep,modeltime(1:3),
     |      timing_type,secs_per_step
        else
          write(6,"('Step ',i6,' of ',i6,' mtime=',3i3)")
     |      istep,nstep,modeltime(1:3)
        endif
      endif ! time2print
!
! Update calendar day and sun's location and declination:
! (sub advance_day also updates sfeps)
!
      if (calendar_advance > 0) call advance_day
!     call sunloc(iday,secs) ! locate sun's longitude
!     if (debug) write(6,"('advance after sunloc')")
!
! Update lower boundaries of t,u,v,z.
      call tuvz_lbc
!
! Calculate vc, barm, and z:
      do j=lat0,lat1
        call addfld('Z_adv0',' ','km',
     |    z(:,lon0:lon1,j,itp)*1.e-5,'lev',1,nlevp1,
     |    'lon',lon0,lon1,j)
      enddo
      z(:,:,:,itp) = 0.
      call addiag(
     |  tn   (levd0,lond0,latd0,itp),   ! in
     |  o1   (levd0,lond0,latd0,itp),   ! in
     |  co   (levd0,lond0,latd0,itp),   ! in
     |  n2   (levd0,lond0,latd0,itp),   ! in
     |  vn   (levd0,lond0,latd0,itp),   ! in
     |  vc   (levd0,lond0,latd0,itp),   ! out
     |  barm (levd0,lond0,latd0),       ! out
     |  co2  (levd0,lond0,latd0),       ! out (3d array in fields.F)
     |  z    (levd0,lond0,latd0,itp),   ! out (note itp is output)
     |  lon0,lon1,1,nlevp1,lat0,lat1)
      if (debug) write(6,"('advance after addiag')")
!
! Define boundary latitudes -1,0 across the south pole, and
!   nlat+1,nlat+2 across the north pole:
!
#ifdef MPI
      call mp_polelats(itp)
      if (debug) write(6,"('advance after mp_polelats')")
!
! Update ghost cell boundary lats and lons.
! Bndlats is called first to exchange lat0-1,2 and lat1+1,2 at lon0->lon1,
!   then bndlons is called  to exchange lon0-1,2 and lon1+1,2 at lat0-2->lat1+2).
!
      f4d(:)%mpi = .true.
      call mp_bndlats(f4d,nf4d,itp)
      call mp_bndlons(f4d,nf4d,itp)
      if (debug) write(6,"('advance after mp_bndlats,mp_bndlons')")
#else
      call mk_polelat(0     ,1     ,itp)
      call mk_polelat(-1    ,2     ,itp)
      call mk_polelat(lat1+1,lat1  ,itp)
      call mk_polelat(lat1+2,lat1-1,itp)
#endif
      do j=lat0,lat1
        call addfld('Z_adv1',' ','cm',
     |    z(:,lon0:lon1,j,itp)*1.e-5,'lev',1,nlevp1,
     |    'lon',lon0,lon1,j)
      enddo
!
! For Z, itc==itp (z at itp was set by addiag):
      z(:,:,:,itc) = z(:,:,:,itp)
      if (debug) write(6,"('advance after z(itc)=z(itp)')")

      do j=lat0,lat1
        call addfld('Z_adv2',' ','cm',
     |    z(:,lon0:lon1,j,itc)*1.e-5,'lev',1,nlevp1,
     |    'lon',lon0,lon1,j)
      enddo
!
! Horizontal diffusion, first step:
      if (debug) write(6,"('advance call hdif1')")
      call hdif1(
     |  tn_nm(levd0,lond0,latd0,itp),
     |  un_nm(levd0,lond0,latd0,itp),
     |  vn_nm(levd0,lond0,latd0,itp),
     |  barm (levd0,lond0,latd0),
     |  1,nlevp1,lon0,lon1,lat0,lat1)
!     if (debug) write(6,"('advance after hdif1')")
!
! Horizontal diffusion, second step:
! hdif2 saves 3d diagnostics f3d (fkldt,u,v,o2,o) at j+1 (0->37)
! hdif2: WRITE (fkldt,u,v,o2,o1) <- READ (tn_nm,un_nm,vn_nm,o2_nm,o1_nm)
! The kld terms will be used later in hdif3.
!
      call hdif2(
     |  tn_nm(levd0,lond0,latd0,itp), ! 4d input
     |  un_nm(levd0,lond0,latd0,itp), ! 4d input
     |  vn_nm(levd0,lond0,latd0,itp), ! 4d input
     |  o1_nm(levd0,lond0,latd0,itp), ! 4d input
     |  co_nm(levd0,lond0,latd0,itp), ! 4d input
     |  kldt ,      ! 3d output
     |  kldu ,      ! 3d output
     |  kldv ,      ! 3d output
     |  kldo1,      ! 3d output
     |  kldco,      ! 3d output
     |  1,nlevp1,lon0,lon1,lat0,lat1)
!     if (debug) write(6,"('advance after hdif2')")
!
! Call main driver (formerly sub dynamics)
      call driver(nstep,istep)

!     do j=lat0,lat1
!       call addfld('Z_adv2',' ','cm',
!    |    z(:,lon0:lon1,j,itc)*1.e-5,'lev',1,nlevp1,
!    |    'lon',lon0,lon1,j)
!     enddo
!
!  -------------------------------------------------------------------------
! 9/06 btf Temporary during vtgcm development:
! 4/08 swb Updated for comp_o2, comp_n4s, comp_n2d
! Define 4-d field at itc from itp, unless the field has been updated 
!   at itc by calls in the driver. As fields are updated, they are added 
!   to this conditional, and the assignment of itc from itp is not executed, 
!   i.e., the updated field will be used in the next timestep.
!  -------------------------------------------------------------------------
!
      do i=1,nf4d
        if (trim(f4d(i)%short_name) /= 'Z'    .and.  ! addiag.F
     |      trim(f4d(i)%short_name) /= 'N2'   .and.  ! compn2.F
     |      trim(f4d(i)%short_name) /= 'OMEGA'.and.  ! swdot.F (W)
     |      trim(f4d(i)%short_name) /= 'TN'   .and.  ! dt.F
     |      trim(f4d(i)%short_name) /= 'TN_NM'.and.  ! dt.F
     |      trim(f4d(i)%short_name) /= 'UN'   .and.  ! duv.F
     |      trim(f4d(i)%short_name) /= 'UN_NM'.and.  ! duv.F
     |      trim(f4d(i)%short_name) /= 'VN'   .and.  ! duv.F
     |      trim(f4d(i)%short_name) /= 'VN_NM'.and.  ! duv.F
     |      trim(f4d(i)%short_name) /= 'O1   '.and.  ! comp.F
     |      trim(f4d(i)%short_name) /= 'O1_NM'.and.  ! comp.F
     |      trim(f4d(i)%short_name) /= 'CO   '.and.  ! comp.F
     |      trim(f4d(i)%short_name) /= 'CO_NM'.and.  ! comp.F
     |      trim(f4d(i)%short_name) /= 'O2  ' .and.  ! comp_o2.F
     |      trim(f4d(i)%short_name) /= 'O2_NM'.and.  ! comp_o2.F
     |      trim(f4d(i)%short_name) /= 'N4S  '.and.  ! comp_n4s.F
     |      trim(f4d(i)%short_name) /= 'N4S_NM'.and. ! comp_n4s.F
     |      trim(f4d(i)%short_name) /= 'NO  ' .and.  ! comp_no.F
     |      trim(f4d(i)%short_name) /= 'NO_NM'.and.  ! comp_no.F
     |      trim(f4d(i)%short_name) /= 'SO2 ' .and.  ! comp_so2.F
     |      trim(f4d(i)%short_name) /= 'SO2_NM'.and. ! comp_so2.F
     |      trim(f4d(i)%short_name) /= 'SO  ' .and.  ! comp_so.F
     |      trim(f4d(i)%short_name) /= 'SO_NM'.and.  ! comp_so.F
     |      trim(f4d(i)%short_name) /= 'N2D'  .and.  ! comp_n2d.F
     |      trim(f4d(i)%short_name) /= 'O2P'  .and.  ! ionden.F
     |      trim(f4d(i)%short_name) /= 'OP'   .and.  ! ionden.F
     |      trim(f4d(i)%short_name) /= 'NE'          ! ionden.F
     |     ) then
          f4d(i)%data(:,:,:,itc) = f4d(i)%data(:,:,:,itp)
        else
!         write(6,"('advance: istep=',i4,': using updated ',a)") 
!    |      istep,trim(f4d(i)%short_name)
        endif
      enddo
!
! Rotate f4d fields with time (sub rotf3d is in rotfld.F):
! subroutine rotf3d(f,lev0,lev1,lon0,lon1,lat0,lat1,step,istep)
!
!     do i=1,nf4d
!       call rotf3d(f4d(i)%data(:,lon0:lon1,lat0:lat1,itc),
!    |    f4d(i)%short_name,1,nlevp1,lon0,lon1,lat0,lat1,step,
!    |    istep,i==nf4d)
!     enddo
!
! Save OMEGA to secondary history:
!     ix = 0
!     do i=1,nf4d
!       if (trim(f4d(i)%short_name)=='OMEGA') ix=i
!     enddo
!     if (ix > 0) then
!       do j=lat0,lat1
!         call addfld('W_f4d','W_f4d',' ',
!    |      f4d(ix)%data(:,lon0:lon1,j,itc),'lev',1,nlevp1,
!    |      'lon',lon0,lon1,j)
!       enddo
!     else
!       write(6,"('>>> advance: could not find f4d index to OMEGA')")
!     endif

!
! Write output history if necessary (output writes data%(:,:,:,itc)).
!
      time2write = wrhist(istep,modeltime,
     |  wrprim, save_prim, newseries_prim, iseries_prim,
     |  wrsech, save_sech, newseries_sech, iseries_sech)
!
! History i/o timing is not included in STEP segment:
      call timer(time0_step,time1,'STEP',1,0) ! suspend step timing for i/o
      time1_step = time1
!
      time1_phist = 0.
      time1_shist = 0.
!
! If its time to write a history, root task must gather subdomain
! data from slave tasks before writing the history. This is done
! by sub mp_gather2root. Only root task writes to the history.
!
#ifdef MPI
      if (time2write.and.ntask > 1) then
        if (wrprim) call timer(time0_phist,time1_phist,'PHIST',0,0)
        call mp_gather2root(itc,'prim')
        if (debug) write(6,"('advance after mp_gather2root')")
!
! Comment this mp_gather2root call if sech fields were written by
!   the master task only (e.g. from serial dynamo).
! This call can also be commented out if only progostics are saved
!   on secondary histories (i.e., addfsech was not called)
!
        if (wrsech) then
          call timer(time0_shist,time1_shist,'SHIST',0,0)
          call mp_gather2root(itc,'sech')
        endif
      endif
#else
      call timer(time0_phist,time1_phist,'PHIST',0,0) ! non-MPI
      call timer(time0_shist,time1_shist,'SHIST',0,0) ! non-MPI
!
! Update foutput if serial non-MPI run:
!
      if (time2write) then
        do i=1,nf4d
          foutput(:,lon0:lon1,lat0:lat1,i) =
     |      f4d(i)%data(:,lon0:lon1,lat0:lat1,itc)
        enddo
        tlbc_glb(lon0:lon1,lat0:lat1) = tlbc(lon0:lon1,lat0:lat1)
        ulbc_glb(lon0:lon1,lat0:lat1) = ulbc(lon0:lon1,lat0:lat1)
        vlbc_glb(lon0:lon1,lat0:lat1) = vlbc(lon0:lon1,lat0:lat1)
      endif
#endif
!
! Write the history to output file:
      if (mytid==0) then
        call outhist(istep,modeltime)
        if (debug) write(6,"('advance after outhist')")
      endif
#ifdef MPI
      call timer(time0_step,time1,'STEP',0,0) ! resume step timing
      if (wrprim) call timer(time0_phist,time1_phist,'PHIST',1,0)
      if (wrsech) call timer(time0_shist,time1_shist,'SHIST',1,0)
      if (timing%level >= 2.and.time2print(nstep,istep)) then
        write(6,"('Time in PHIST =        ',f12.2,
     |    ' Advance: step ',i5)") time1_phist,istep
        write(6,"('Time in SHIST =        ',f12.2,
     |    ' Advance: step ',i5)") time1_shist,istep
      endif
!
! Each mpi task must receive the 2 lats before its first
! updated lat (lat0-1,lat0-2), and the 2 lats after its last
! updated lat (lat1+1,lat2+2). Each task must also send its
! first 2 (lat0,lat0+1) and last 2 (lat1,lat1-1) updated lats
! to the appropriate "adjacent" tasks.
!
      call timer(time0_prep,time1,'PREP',0,0) ! resume prep timing
      f4d(:)%mpi = .true.
      call mp_bndlats(f4d,nf4d,itc)
      if (debug) write(6,"('advance after mp_bndlats')")
      call mp_bndlons(f4d,nf4d,itc)
      if (debug) write(6,"('advance after mp_bndlons')")
!
! Periodic points for all updated fields:
      call mp_periodic_f4d(itc)
#else
! Non-mpi serial run:
      call timer(time0_step,time1,'STEP',0,0) ! resume step timing (non-MPI)
      call mk_polelat(0     ,1     ,itc)
      call mk_polelat(-1    ,2     ,itc)
      call mk_polelat(lat1+1,lat1  ,itc)
      call mk_polelat(lat1+2,lat1-1,itc)
      call set_periodic_f4d(itc)
#endif
!
! Swap field data time indices, so current updated data becomes previous
! step data for next step:
      itmp = itp
      itp = itc
      itc = itmp
!
      call timer(time0_prep,time1,'PREP',1,0) ! end prep timing
      time1_prep = time1_prep+time1
      call timer(time0_step,time1,'STEP',1,0) ! end step timing
      time1_step = time1_step+time1
      if (timing%level >= 2.and.time2print(nstep,istep)) then
        write(6,"('Time in PREP =         ',
     |    f12.2,' Advance: step ',i5)") time1_prep,istep
        write(6,"('Time in STEP =         ',
     |    f12.2,' Advance: step ',i5)") time1_step,istep
      endif
      secs_per_step = time1_step
!
! Return for next time step:
      if (istep < nstep) then
#ifdef MPI
        if (ntask > 1) then
          call mpi_barrier(MPI_COMM_WORLD,ier)
          if (debug) write(6,"('advance after end of step barrier')")
        endif
#endif
        if (debug) write(6,"('advance end iteration for step ',i4)")
     |    istep
        goto 100
      endif
!
#ifdef VT
!     code = 115 ; state = 'advance' ; activity='ModelCode'
      call vtend(115,ier)
#endif
!
      end subroutine advance
!-----------------------------------------------------------------------
      subroutine advance_day
!
! Advance calendar day if needed. Also update sfeps.
!
      use init_module,only: iter,iyear,iday,sfeps,sundec,sin_sundec,
     |  cos_sundec
      use cons_module,only: dt,pi
      implicit none
!
! Local:
      integer :: idayit,idayprev,iyearprev,iyr4,iyr100,lpyr,ienda
      real :: delta,theta0
!
      idayit = iter*int(dt)/86400
      if (idayit*86400 == iter*int(dt)) then
        idayprev = iday
        iyearprev = iyear
        iday = iday + 1
!
! lpyr = 1(0) if is (not) a leap year
        iyr4 = iyear/4
        iyr100 = iyear/100
        lpyr = 0
        if (iyr4*4 == iyear .and. iyr100*100 /= iyear) lpyr=1
        ienda = 365 + lpyr
        if (iday > ienda) then
          iyear = iyear + 1
          iday = iday - ienda
        endif                          !   for past year's end
!
! Recalculate sun's declination
!       delta = atan(tan(23.5*pi/180.)*sin(2.*pi*real(iday-80)/365.))
        delta = 0.0
        sin_sundec = sin(delta) ! C(95)
        cos_sundec = cos(delta) ! C(96)
!
! Update sfeps:
! sfeps is 6% variation in solar output over a year
!   caused by the orbital eccentricity.
!
        theta0 = 2.*pi*real(iday)/365.
        sfeps = 1.000110+0.034221*cos(theta0)+0.001280*sin(theta0)
     1        +0.000719*cos(2.*theta0)+0.000077*sin(2.*theta0)
!
        write(6,"('Advancing day (previous,present)=',4i5,' sfeps=',
     |   e12.4)") idayprev,iyearprev,iday,iyear,sfeps
      endif
      end subroutine advance_day
!-----------------------------------------------------------------------
!     logical function time2print(nstep,istep)
!     implicit none
!       integer,intent(in) :: nstep,istep
!       time2print = .false.
!       if (nstep <= 100 .or. (nstep > 100 .and. mod(istep,10)==0))
!    |    time2print = .true.
!       time2print = .true.   ! print every step
!     end function time2print
!-----------------------------------------------------------------------

