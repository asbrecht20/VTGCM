! 
      subroutine settei(tn,z,te,ti,lev0,lev1,lon0,lon1,lat)
!
! ---------------------------------------------------------------------
! ** Calculate electron temperatures.
! Use Theis and Brace (1993) Te formulation for PVO entry (nightside)
!     and modify it slightly to reproduce Fox & Sung(2001) Te.
!  (a)  At/Below dayside ionospheric peak, Te = Tn is fine
!  (b)  Over 140-165 km use Theis and Brace scheme modified by -0.27
!  (c)  Above 165 km use Theis and Brace scheme modified by -0.061
! Use Fox & Sung(2001) Te formulation and test impacts upon DR reactions.
!  (a)  At/Below dayside ionospheric peak, Te = Tn is fine
!  (b)  Over 140-165 km use formulation #1
!  (c)  Above 165 km use formulation #2
! Use Theis and Brace (1993) Te formulation for PVO entry (nightside)
!     and modify to activate above 160 km (Like old VTGCM)
!  (a)  Over <160 km  Te = Tn
!  (b)  Above 160 km use Theis and Brace scheme 
! ** Calculate ion temperatures.
! Assume that Ti=Tn at all altitudes for now.
! ---------------------------------------------------------------------
!
      use params_module,only: nlevp1,nlonp4
      use addfld_module,only: addfld
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,      ! neutral temperature (deg K)
     |   z       ! absolute height (cm) from addiag.F
!
! Output args:
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) ::
     |  te,      ! output electron temperature (deg K) 
     |  ti       ! output ion temperature (deg K)
!
! Local:
      integer :: k,i
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  zkm,    ! absolute height (in km) 
     |  logte   ! log10te (log of Te)
!
! ---------------------------------------------------------
! Set ion temperature output. 
! Update later from Fox and Sung (2001)
! -- Ti = Tn assumption is fine below 160 km.
! -- Ti = Tn will adversely impact N2+ and NO+ S and L
!         at the highest altitudes (<=20% level).
! ---------------------------------------------------------
          ti(:,:) = tn(:,:) 
!
! ---------------------------------------------------------
!  zkm derived from geopotential height (from addiag.F)
! ---------------------------------------------------------
      do i=lon0,lon1
        do k=lev0,lev1-1
          zkm(k,i) = z(k,i)/1.0E+05
        enddo
      enddo
!  -------------------------------------------------------
!  Te formulation from Theis and Brace (1993)
!   GRL, 20, 2719, (1993). PVO Entry (niteside)
!  -------------------------------------------------------
      do k=lev0,lev1-1
        do i=lon0,lon1
! ----------------------------------------------------------
!  ** Te formulation from Theis and Brace (1993); Modified
!  ** Fox and Sung (2001) modification
!  Te = Tn below 140 km (at/below ionospheric peak)
!       if (zkm(k,i) <= 140.) then
!           logte(k,i)= log10(tn(k,i))
!        elseif (zkm(k,i) > 140. .and. zkm(k,i) <= 165.) then
!  Te = formulation over 140-165 km
!          logte(k,i)= 3.201 - 1921.9/(zkm(k,i)-98.078)**2.+
!    |                 8.5257*zkm(k,i)*1.0E-04
!        else
!  Te = formulation above 165 km
!          logte(k,i)= 3.410 - 1921.9/(zkm(k,i)-98.078)**2.+
!    |                 8.5257*zkm(k,i)*1.0E-04
!        endif
!       enddo !  i=lon0,lon1
! ----------------------------------------------------------
!  ** Theis and Brace (1993) formulation (straight)
!  Te = Tn below 140 km (at/below ionospheric peak)
 	 if (zkm(k,i) <= 140.) then
 	   logte(k,i)= log10(tn(k,i))
         else
!  Te = formulation above 140 km
           logte(k,i)= 3.471 - 1921.9/(zkm(k,i)-98.078)**2.+
     |                 8.5257*zkm(k,i)*1.0E-04
         endif
! --------------------------------------------------------
!  ** Te formulation from Theis and Brace (1993); Modified
!   GRL, 20, 2719, (1993). PVO Entry (niteside)
! --------------------------------------------------------
!  Te = Tn below 160 km (at/below/above ionospheric peak)
!	 if (zkm(k,i) <= 160.) then
!	   logte(k,i)= log10(tn(k,i))
!        else
!  Te = formulation above 160 km
!         logte(k,i)= 3.471 - 1921.9/(zkm(k,i)-98.078)**2.+
!    |                 8.5257*zkm(k,i)*1.0E-04
!        endif
! ----------------------------------------------------------
        enddo !  i=lon0,lon1
      enddo !  k=lev0,lev1-1 
! ---------------------------------------------------------
! Te extracted from log10(Te)
! ---------------------------------------------------------
      do i=lon0,lon1
        do k=lev0,lev1-1
          te(k,i) = 10.**logte(k,i) 
        enddo
      enddo
! ---------------------------------------------------------
      call addfld('TE',' ',' ',te(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('TI',' ',' ',ti(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
! ---------------------------------------------------------
      end subroutine settei
