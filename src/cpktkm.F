 
      subroutine cpktkm(tn,o1,co,co2,n2,fcp,fkt,fkm,
     |  lev0,lev1,lon0,lon1,lat)
!
! Define diagnostics CP, KT, and KM.
!
      use cons_module,only: gask,rmassinv_o1,rmassinv_co,
     |  rmassinv_n2,rmassinv_co2
      use addfld_module,only: addfld
      implicit none
!
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,  ! neutral temperature (deg K)
     |  o1,
     |  co,
     |  co2, ! 1-o-co-n2 from addiag.F
     |  n2
! Output args:
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) ::
     |  fcp, ! specific heat at constant pressure (ergs/deg/gm)
     |  fkt, ! molecular diffusion (ergs/cm/deg/sec)
     |  fkm  ! molecular viscosity (gm/cm/sec)
!
! Local:
      real :: 
     |  cmrf(8) = (/135.8,185.7,230.4,271.0,308.3,343.3,373.7,406.1/),
     |  com (8) = (/0.3966,0.7692,1.0776,1.340,1.574,1.787,1.986,2.172/)
      real :: ttot,cokm,prco,co2kt,cokt,rrco2,tt,cpco2,co2km,crn
      integer :: k,i,is
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  mbar,    ! mean molecular weight (s1)
     |  po,pco,pco2,pn2 ! mixing ratios (s2,s3,s4,s5)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          mbar(k,i) = 1./(o1(k,i)*rmassinv_o1+co(k,i)*rmassinv_co+
     |      n2(k,i)*rmassinv_n2+co2(k,i)*rmassinv_co2)

          po(k,i)  = mbar(k,i)*o1(k,i)*rmassinv_o1
          pco(k,i) = mbar(k,i)*co(k,i)*rmassinv_co
          pn2(k,i) = mbar(k,i)*n2(k,i)*rmassinv_n2
          pco2(k,i)= mbar(k,i)*co2(k,i)*rmassinv_co2

          fcp(k,i) = gask*(pco2(k,i)*(3.5*rmassinv_co2)+
     |      (pn2(k,i)+pco(k,i))*(3.5*rmassinv_co)+
     |      po(k,i)*(2.5*rmassinv_o1))

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C---- KT, KM FORMULATION TAKEN FROM VTGCM INPUT DATASET CODE
C     ( BANKS AND KOCKARTS )
C---- KM=((PO*3.9)+(PN2*3.42))*TT*1.E-06 +(PCO*COKM)+(PCO2*CO2KM)
C---- KT=((PO*75.9)+(PN2*56.))*TT +(PCO*COKT)+(PCO2*CO2KT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          ttot = tn(k,i)
          tt = ttot**0.69
!
! co2 factors:
          is = int((ttot-173.3)/100.)
          if (is <= 1) is = 1
          if (is >= 7) is = 7
          rrco2 = gask*rmassinv_co2
          if (ttot < 500.) crn = 1.64-(ttot-500.)*2.5e-4
          co2km=cmrf(is)+(cmrf(is+1)-cmrf(is))*(ttot-(is*100.+73.3))*
     |      0.01
          co2km=co2km*1.e-06
          cpco2=3.5*gask*rmassinv_co2
          co2kt=(cpco2-rrco2)*co2km*crn
!
! co factors:
          is=int(ttot/100.)
          if (is <= 1) is=1
          if (is >= 7) is=7
          if (ttot > 400.) prco=0.72
          if (ttot > 300. .and. ttot < 400.) prco = 0.73-(ttot-350.)*
     |      1.5e-04
          if (ttot < 300.) prco=0.75-(ttot-250.)*2.6e-04

          cokm=com(is)+(com(is+1)-com(is))*(ttot-100.*is)*0.01
          cokm=cokm*1.65e-04
          cokt=(3.5*gask*cokm)/(28.*prco)
!
! Total mixture kt and km formulation
! po,pco,pco2,pn2 ! mixing ratios (s2,s3,s4,s5)
          fkm(k,i) = (pn2(k,i)*3.42+po(k,i)*3.9)*1.0e-06*tt+
     |      pco2(k,i)*co2km+pco(k,i)*cokm
          fkt(k,i) = (po(k,i)*75.9+pn2(k,i)*56.)*tt+pco(k,i)*cokt+
     |      pco2(k,i)*co2kt

        enddo ! i=lon0,lon1
      enddo ! k=lev0,lev1-1
!
      do i=lon0,lon1
        do k=lev1,lev0+1,-1
          fcp(k,i) = .5*(fcp(k,i)+fcp(k-1,i))
          fkm(k,i) = .5*(fkm(k,i)+fkm(k-1,i))
          fkt(k,i) = .5*(fkt(k,i)+fkt(k-1,i))
        enddo
      enddo
      do i=lon0,lon1
        fcp(lev0,i) = 2.*fcp(lev0,i)-fcp(lev0+1,i)
        fkm(lev0,i) = 2.*fkm(lev0,i)-fkm(lev0+1,i)
        fkt(lev0,i) = 2.*fkt(lev0,i)-fkt(lev0+1,i)
      enddo

!     call addfld('CP',' ',' ',fcp(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('KM',' ',' ',fkm(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('KT',' ',' ',fkt(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)

      end subroutine cpktkm
