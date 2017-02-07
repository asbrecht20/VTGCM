 
      subroutine compn2(tn,mbar,psn2,lev0,lev1,lon0,lon1,lat)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ****     N2 COMPOSITION CALC FOR NEW DTUV BACKGRD :  NN2 SLOT
C     ****     DIAGNOSTIC ONLY. MUST PROCEED COMP CALCULATION.
C     ****     TOTAL CALCULATION ON MID-LEVELS!
C     ****     ONCE EACH NSTEP
C     ****     --FILTER2 (ON) TO SMOOTH LONGITUDE WAVES IN N2 (11/04/99)  OK!
C     ****     --NJNP SLOT USED FOR UPDATED N2 (10/28/99)  OK!
C     ****     Latest Version:  10/26/06 Bougher     venus xmixn2 used
C     ****     Latest Version:  05/14/06 Bougher     izpt = 31 (z = -1.0 homopause)
C     ****     Latest Version:  05/14/06 Bougher     xmixn2 = 0.043
C     ****     Latest Version:  05/15/06 Bougher     xmixn2 = 0.069
C     ****     Latest Version:  05/30/06 Bougher     izpt = 33 (z = 0.0 homopause)
C     ****     Latest Version:  05/30/06 Bougher     xmixn2 = 0.100
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      use params_module,only: dz,zmbot,nlevp1
      use cons_module,only: kut,rmass,p0,gask,grav,boltz,rmass_n2
      use addfld_module,only: addfld
      use chapman_module,only: idn_lonlat
      implicit none
!
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in)::
     |  tn,   ! neutral temperature
     |  mbar  ! mean molecular weight (mbar from addiag)
!
! Output args:
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) ::
     |  psn2  ! n2 (mmr)
!
! Local:
!     real,parameter :: xmixn2=0.041
!     real,parameter :: xmixn2=0.043
!     real,parameter :: xmixn2=0.069
      real,parameter :: xmixn2=0.100
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat)::
     |  n2ref                                ! Ref. mmr for testing
      real,dimension(nlevp1)  :: n2refd,  ! Ref. vmr for testing
     |  n2refn,mbarrefd,mbarrefn
      integer :: i,k,izpt
      real :: szp
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  mbarm,  ! mbar at midpoints        (s5)
     |  nmbar,  ! nmbar                    (s1)
     |  scht,   ! mean scale height        (s6)
     |  n2scht, ! n2 pressure scale height (s7) 
     |  xn2     ! local n2                 (s3)

      if (dz == 0.5) then
        izpt=33
      elseif (dz == 0.25) then
        izpt=66
      else
        write(6,"('>>> compn2: Unknown dz=',f10.3)") dz
        call shutdown('compn2: unknown dz')
      endif


! Day

!      n2refd(:) = (/
!    | 2.2915e-05,3.4037e-06,2.6878e-06,4.9127e-06,9.2205e-06,
!    | 1.8368e-05,3.7627e-05,7.6915e-05,1.5689e-04,3.3004e-04,
!    | 6.9942e-04,1.4794e-03,3.1012e-03,6.4967e-03,1.3549e-02,
!    | 2.8098e-02,4.0600e-02,4.1940e-02,4.3029e-02,4.4301e-02,
!    | 4.5518e-02,4.7109e-02,4.8950e-02,5.0799e-02,5.3007e-02,
!    | 5.5256e-02,5.7925e-02,6.0959e-02,6.3953e-02,6.7018e-02,
!    | 6.9955e-02,7.2651e-02,7.5051e-02,7.6837e-02,7.7200e-02,
!    | 7.6293e-02,7.3672e-02,6.9543e-02,6.4174e-02,5.7274e-02,
!    | 4.9395e-02,4.0716e-02,3.2203e-02,2.4688e-02,
!    | 2.4688e-02,2.4688e-02,2.4688e-02,2.4688e-02,2.4688e-02,
!    | 2.4688e-02,2.4688e-02,2.4688e-02,2.4688e-02,2.4688e-02,
!    | 2.4688e-02,2.4688e-02,2.4688e-02,2.4688e-02,2.4688e-02,
!    | 2.4688e-02,2.4688e-02,2.4688e-02,2.4688e-02,2.4688e-02,
!    | 2.4688e-02,2.4688e-02,2.4688e-02,2.4688e-02,2.4688e-02/)
!
!      mbarrefd(:) = (/
!    | 44.,44.,44.,44.,44.,44.,44.,44.,44.,44.,
!    | 43.9,43.9,43.9,43.9,43.8,43.5,43.3,43.3,
!    | 43.3,43.3,43.2,43.2,43.2,43.1,43.0,42.9,
!    | 42.7,42.5,42.2,41.8,41.2,40.3,39.2,37.7,
!    | 35.9,33.7,31.4,29.0,26.6,24.4,22.4,20.6,
!    | 20.6,20.6,20.6,20.6,20.6,20.6,20.6,20.6,
!    | 20.6,20.6,20.6,20.6,20.6,20.6,20.6,20.6,
!    | 20.6,20.6,20.6,20.6,20.6,20.6,20.6,20.6,
!    | 20.6,20.6,20.6/)
!
! Night
!      n2refn(:) = (/
!    | 2.2915e-05,3.4037e-06,2.6878e-06,4.9127e-06,9.2205e-06,
!    | 1.8368e-05,3.7627e-05,7.6915e-05,1.5689e-04,3.3563e-04,
!    | 7.3580e-04,1.6155e-04,3.5764e-03,7.8137e-03,1.7171e-02,
!    | 3.4163e-02,4.0061e-02,4.0822e-02,4.1686e-02,4.2758e-02,
!    | 4.3938e-02,4.5216e-02,4.6687e-02,4.8311e-02,5.0107e-02,
!    | 5.2345e-02,5.4829e-02,5.7500e-02,6.0551e-02,6.3930e-02,
!    | 6.7457e-02,7.1440e-02,7.5949e-02,8.0693e-02,8.5386e-02,
!    | 8.9875e-02,9.4445e-02,9.8040e-02,9.7861e-02,9.6880e-02,
!    | 9.2571e-02,8.2073e-02,7.1298e-02,5.5677e-02,4.2086e-02,
!    | 2.9856e-02,2.0642e-02,1.3961e-02,9.2085e-03,5.9392e-03,
!    | 3.6518e-03,2.1612e-03,1.1744e-03,5.6280e-04,2.2196e-04,
!    | 2.2196e-04,2.2196e-04,2.2196e-04,2.2196e-04,2.2196e-04,
!    | 2.2196e-04,2.2196e-04,2.2196e-04,2.2196e-04,2.2196e-04,
!    | 2.2196e-04,2.2196e-04,2.2196e-04,2.2196e-04/)
!
!      mbarrefn(:) = (/
!    | 44.,44.,44.,44.,44.,44.,44.,44.,44.,44.,
!    | 43.9,43.9,43.8,43.7,43.6,43.3,43.1,43.2,
!    | 43.2,43.2,43.1,43.1,43.2,43.2,43.0,42.7,
!    | 42.6,42.4,42.2,41.8,41.4,40.8,40.0,38.5,
!    | 37.1,34.5,32.0,29.7,27.7,26.2,25.7,26.6,
!    | 26.6,26.6,26.6,26.6,26.6,26.6,26.6,26.6,
!    | 26.6,26.6,26.6,26.6,26.6,26.6,26.6,26.6,
!    | 26.6,26.6,26.6,26.6,26.6,26.6,26.6,26.6,
!    | 26.6,26.6,26.6/)


!** convert vmr to mmr
!       do i=lon0,lon1   !1,nlonp4
!         do k=lev0,lev1   !1,nlevp1
!    Day-Night Variations Specified
!           if (idn_lonlat(i,lat)==1) then
!             n2ref(k,i,lat) = n2refd(k)*rmass_n2/mbarrefd(k)
!           else
!             n2ref(k,i,lat) = n2refn(k)*rmass_n2/mbarrefn(k)
!           endif
!         enddo
!       enddo
!
! Save n2 that was read from data file by read_n2file:
!     call addfld('N2_DAT',' ',' ',n2dat(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('COMPN2_T',' ',' ',tn(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! CALC ON MID-LEVELS THROUGHOUT!
!
! mbarm = mbar at midpoints:
      do i=lon0,lon1
        do k=lev0,lev1-1
          mbarm(k,i) = 0.5*(mbar(k,i)+mbar(k+1,i))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!     call addfld('MBARM',' ',' ',mbarm(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! nmbar = n*mbar
      do i=lon0,lon1
        do k=lev0,lev1-1
          szp = zmbot + (0.5*dz) + (k-1)*dz
          nmbar(k,i) = p0*exp(-szp)/(boltz*tn(k,i))*mbarm(k,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!     call addfld('NMBAR',' ',' ',nmbar(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! scht = mean scale height
! n2scht = n2 pressure scale height 
      do i=lon0,lon1
        do k=lev0,lev1-1
          scht(k,i) = gask*tn(k,i)/(mbarm(k,i)*grav)
          n2scht(k,i) = gask*tn(k,i)/(rmass_n2*grav)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!     call addfld('SCHT',' ',' ',scht(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('N2SCHT',' ',' ',n2scht(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! xn2 = n2 at and below z=zpt
      do i=lon0,lon1
        do k=1,izpt
          xn2(k,i) = nmbar(k,i)/mbarm(k,i)*xmixn2
        enddo
        do k=izpt+1,lev1-1
          xn2(k,i) = xn2(k-1,i)*exp(-(dz*scht(k,i)/n2scht(k,i)))
        enddo
      enddo
!     call addfld('XN2',' ',' ',xn2(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! psi n2 output:
      do i=lon0,lon1
        do k=lev0,lev1-1
          psn2(k,i) = xn2(k,i)*rmass_n2/nmbar(k,i)
        enddo ! k=lev0,lev1-1
        psn2(lev1,i) = psn2(lev1-1,i)
      enddo ! i=lon0,lon1
! Overwrite with Ref Atm to test cool 
!     do i=lon0,lon1
!       do k=lev0,lev1-1
!         psn2(k,i) = n2ref(k,i,lat)
!       enddo ! k=lev0,lev1-1
!     enddo ! i=lon0,lon1
!     call addfld('PSN2',' ',' ',psn2(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Filter_sp is called by driver.F, after this routine returns.
!
      end subroutine compn2
!-----------------------------------------------------------------------
