!
      program tgcm
!
! Main program for tgcm models.
! -- latest revision : btf (8/1/07): call qirgrid_init
! -- latest revision : swb (9/25/07): call coolref_init
! -- latest revision : asb (5/19/08): call pepdgrid_init
!
      use input_module,only: input,dynamo,step,datavol,pid,hpss_path
      use dispose_module,only: init_dispose
      use init_module,only: init,iyear
      use hist_module,only: nstep
      use lbc,only: init_lbc
      use timing_module,only: timer,timer_report,init_timer
      use qrj_module,only: init_qrj,init_euvac,init_sigmas,init_sflux
      use qirgrid_init_module,only: qirgrid_init
      use coolref_init_module,only: coolref_init
      use pepdgrid_init_module,only: pepdgrid_init
#ifdef MPI
      use mpi_module,only: mp_init,mp_close,mytid,ntask,mp_distribute
#endif
      implicit none
!
! Local:
      integer :: ier,icount_tgcm,icount_apex,nsecs
      real :: elapsed,cpu1,cpu2
      real ::
     |  time0,     time1,
     |  time0_run, time1_run,
     |  time0_init,time1_init
!
! Report starting time and other info to stdout:
      call startup_message
!
! Initialize timing for entire run:
#ifndef IRIX
      call cpu_time(cpu1)
#endif
!
! Initialize timing (must be called after mp_init):
      call init_timer
!
! Get user input. If MPI job, init mpi and set up 2-d decomposition
!   across tasks. Start timing contains an mpi_barrier, so must be
!   called after mp_init.
!
#ifdef MPI
      call mp_init
#endif
      call timer(time0_run,time1_run,'RUN',0,0)    ! start run timing
      call timer(time0_init,time1_init,'INIT',0,0) ! start init timing
      call getpid(pid)
      write(6,"('Process ID (pid) = ',i8)") pid
!
#ifdef MPI
      call input(mytid,ntask) ! get user input (mpi)
      call mp_distribute
#else
      call input(-1,0)    ! get user input (non-mpi)
#endif
      if (mytid==0.and.len_trim(hpss_path) > 0) call init_dispose
!
! Do initialization:
      call init
!
! Initialization for qrj:
! (this is not in init_module to avoid circular dependency between
!  init_module and qrj_module)
      call init_qrj
      call init_euvac
      call init_sigmas
      call init_sflux
      call qirgrid_init
      call coolref_init
      call pepdgrid_init
!
! Read source history:
      call readsource(ier)
!
! Read data volume(off): old DRM86 3-D and 1-D arrays.
!     call read_ncdata(datavol)
!
! Call apex code if doing dynamo.
!     if (dynamo > 0) call apxparm(real(iyear))
!
! Set up magnetic field data (read from magfield file if dynamo<=0)
!     call magfield(dynamo)
!
! Set lower boundary parameters (including Hough mode functions):
      call init_lbc

      call timer(time0_init,time1_init,'INIT',1,0) ! end init timing
!
! Advance the model (timing in main time-step loop is done in advance):
      call advance
!
! Report to stdout:
      write(6,"(' ')")
#ifdef MPI
      write(6,"('MPI run with ntask = ',i3)") ntask
#endif
      write(6,"('nstep=',i5,' step=',i5)") nstep,step
      nsecs = nstep*step
      write(6,"('Model simulation time = ',i8,' secs ',/,
     |  '  (minutes=',f8.2,', hours=',f8.2,', days=',f10.6,')')")
     |  nsecs,float(nsecs)/60.,float(nsecs)/3600.,
     |  float(nsecs)/(24.*3600.)
!
! End timing:
!
#ifndef IRIX
      call cpu_time(cpu2)
      write(6,"('Cpu time for run = ',f10.2)") cpu2-cpu1
#endif
      call timer(time0_run,time1_run,'RUN',1,0) ! end total run timing
!
! Finalize mpi:
!
#ifdef MPI
      call mp_close ! mpi.F
#endif
      call timer_report
!
      call final_message
      write(6,"('NORMAL EXIT')")
      end program tgcm
!-----------------------------------------------------------------------
      subroutine startup_message
      character(len=8) ::
     |  rundate,           ! current local date
     |  runtime            ! current local time
      character(len=16) ::
     |  host,              ! host machine
     |  system,            ! operating system of host (from pre-proc macros)
     |  logname            ! user login name

      call datetime(rundate,runtime) 
      call gethostsname(host)
      call setosys(system)
      logname = ' '
      call getenv('LOGNAME',logname)
      if (len_trim(logname)==0) logname = "unknown"
      write(6,"(/,72('='))")
      write(6,"('Begin execution of TIEGCM at ',a,'  ',a)")
     |  rundate,runtime
      write(6,"('  Host    = ',a)") trim(host)
      write(6,"('  System  = ',a)") trim(system)
      write(6,"('  Logname = ',a)") trim(logname)
      write(6,"  (72('='),/)")
      end subroutine startup_message
!-----------------------------------------------------------------------
      subroutine final_message
!
! Print end-of-execution message to stdout with date and time:
!
      use params_module,only: tgcm_version
      character(len=8) ::
     |  rundate,           ! current local date
     |  runtime            ! current local time
      call datetime(rundate,runtime)
      write(6,"('End execution of ',a,' at ',a,'  ',a)")
     |  trim(tgcm_version),rundate,runtime
      end subroutine final_message

