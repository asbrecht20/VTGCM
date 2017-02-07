      module ionden_module
      implicit none
      contains
!-----------------------------------------------------------------------
      subroutine ionden(tn,barm,o1,co,co2,n4s,n2,
     |  n2p,nop,co2p,o2p,op,ne,lev0,lev1,lon0,lon1,lat)
!
! Calculate production and loss terms for ion species.
! Return ion densities as production/loss
! Add Iteration Capability similar to vtgcm2 (1997)
! -- Latest mods : swb (11/15/06) iteration capability, initial conditions
! -- Latest mods : swb (11/16/06) recoding for neutral densities
! -- Latest mods : swb (12/03/07) cutoff ion and electron densities at 0.01 #/cm3
!                  for code, but fix at 1.0 for plotter. Works with comp_oco.F
!                  for no OSRC or COSRC below z = -12.0 for now.
! -- Latest mods:  swb (4/18/08)  drop 0.01 minimum for ions, retain for Ne.
!
      use params_module,only: spval
      use qrj_module,only: qop,qn2p,qco2p,qo2p, ! qnop not available 
     |  mk_xnmbari
      use fields_module,only: tlbc
      use cons_module,only: rmassinv_o1,rmassinv_co,rmassinv_co2,
     |  rmassinv_n4s,rmassinv_n2
      use addfld_module,only: addfld
      use chemrates_module,only: 
     |  r18,r19,r20,r28,             ! constant (scalar parameters)
     |  r21,r27,r29,a13,a14,a15,a17  ! temperature-dependent (k,i,j)
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
!
! Input args: press vs longitude (2d (k,i) on midpoints):
!
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: 
     |  tn,    ! neutral temperature (deg K)
     |  barm,  ! mean molecular weight
     |  o1,    ! o1 (mmr)
     |  co,    ! co (mmr)
     |  co2,   ! co2 (mmr)
     |  n4s,   ! n4s (mmr)
     |  n2     ! n2 (mmr)
!
! Output args (particles/cm3):
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) :: 
     |  n2p,       ! N2+
     |  nop,       ! NO+
     |  co2p,      ! CO2+
     |  o2p,       ! O2+
     |  op,        ! O+
     |  ne         ! electron density output
!
! Local:
      integer :: k,i,k0,k1,i0,i1,it
      real :: test, neold
!
! Production terms in cm3/s, loss terms in 1/s
! Number densities (#/cm3) created from mmr values in input arg list
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  n2p_prod,  n2p_loss,
     |  nop_prod,  nop_loss,
     |  co2p_prod, co2p_loss,
     |  o2p_prod,  o2p_loss,
     |  op_prod,   op_loss,
     |  xnmbari,  ! N*MBAR at interfaces
     |  xnmbar,   ! N*MBAR at midlevels 
     |  tni,      ! tn at interfaces
     |  xo,       ! o1 (den)
     |  xco,      ! co (den)
     |  xco2,     ! co2 (den)
     |  xn4s,     ! n4s (den)
     |  xn2       ! n2 (den)
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  n2pl,       ! log10(N2+)
     |  nopl,       ! log10(NO+)
     |  co2pl,      ! log10(CO2+)
     |  o2pl,       ! log10(O2+)
     |  opl,        ! log10(O+)
     |  nel         ! log10(electron density)
!
! For addfld calls:
      k0=lev0 ; k1=lev1 ; i0=lon0 ; i1=lon1
!
! Calculate tn at interfaces:
      do i=lon0,lon1
      tni(lev0,i) = tlbc(i,lat)
      	do k=lev0+1,lev1-1
         	tni(k,i) = .5*(tn(k-1,i)+tn(k,i))
      	enddo
      tni(lev1,i) = tn(lev1-1,i) ! nlevp1 <- nlev
      enddo
! --------------------------------------------------------------------
!
! Calculate p0*e(-z)*barm/kT (1/(N*MBAR))at interfaces
!     subroutine mk_xnmbari(tni,barm,xnmbari,lev0,lev1,lon0,lon1,lat)
!
      call mk_xnmbari(tni,barm(:,lon0:lon1),
     |  xnmbari,
     |  lev0,lev1,lon0,lon1,lat)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          xo(k,i)   = xnmbari(k,i)*o1(k,i) *rmassinv_o1
          xco(k,i)  = xnmbari(k,i)*co(k,i) *rmassinv_co
          xco2(k,i) = xnmbari(k,i)*co2(k,i)*rmassinv_co2
          xn4s(k,i) = xnmbari(k,i)*n4s(k,i)*rmassinv_n4s
          xn2(k,i)  = xnmbari(k,i)*n2(k,i) *rmassinv_n2
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!     call addfld('XCO2',' ',' ',xco2,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XCO' ,' ',' ',xco ,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XN2' ,' ',' ',xn2 ,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XO'  ,' ',' ',xo  ,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XN4S',' ',' ',xn4s,'lev',k0,k1,'lon',i0,i1,lat)
!
! --------------------------------------------------------------------
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          ne(k,i) = 5.0E+05 ! Initialization (until available in history)
!
! --------------------------------------------------------------------
! Iteration loop for PCE ions and Ne (vtgcm2 scheme as benchmark for elden)
!
          do it=1,100
!
! N2+  (first)
          n2p_prod(k,i) = qn2p(k,i,lat)
          n2p_loss(k,i) = r21(k,i,lat)*xco2(k,i)+r27(k,i,lat)*xo(k,i)+
     |                    a15(k,i,lat)*ne(k,i) 
          n2p(k,i) = n2p_prod(k,i)/n2p_loss(k,i)
!
! CO2+ (second)
          co2p_prod(k,i)= qco2p(k,i,lat)+r21(k,i,lat)*n2p(k,i)*xco2(k,i)
          co2p_loss(k,i)= (r19+r20)*xo(k,i)+a14(k,i,lat)*ne(k,i)
          co2p(k,i) = co2p_prod(k,i)/co2p_loss(k,i)
!
! O+ (third)
          op_prod(k,i) = qop(k,i,lat)+r19*xo(k,i)*co2p(k,i)
          op_loss(k,i) = r18*xco2(k,i)+r29(k,i,lat)*xn2(k,i)
!         op_loss(k,i) = r18*xco2(k,i)+r29*xn2(k,i)
          op(k,i) = op_prod(k,i)/op_loss(k,i)
!
! O2+ (fourth) 
          o2p_prod(k,i) = qo2p(k,i,lat)+r20*xo(k,i)*co2p(k,i)+
     |      r18*op(k,i)*xco2(k,i)
          o2p_loss(k,i) = a13(k,i,lat)*ne(k,i)+r28*xn4s(k,i)
          o2p(k,i) = o2p_prod(k,i)/o2p_loss(k,i)
!
! NO+ (fifth)
          nop_prod(k,i) = r28*o2p(k,i)*xn4s(k,i)+r27(k,i,lat)*n2p(k,i)*
     |      xo(k,i)+r29(k,i,lat)*op(k,i)*xn2(k,i)
!    |      xo(k,i)+r29*op(k,i)*xn2(k,i)
          nop_loss(k,i) = a17(k,i,lat)*ne(k,i)
          nop(k,i) = nop_prod(k,i)/nop_loss(k,i)
!
! Ne (last) and check for values less than 0.01
! Retain simulated values for all ions
	  ne(k,i) = n2p(k,i)+co2p(k,i)+op(k,i)+o2p(k,i)+nop(k,i)
          if (ne(k,i) < 0.01) ne(k,i) = 0.01
!         if (nop(k,i) < 0.01) nop(k,i) = 0.01
!         if (o2p(k,i) < 0.01) o2p(k,i) = 0.01
!         if (op(k,i) < 0.01) op(k,i) = 0.01
!         if (co2p(k,i) < 0.01) co2p(k,i) = 0.01
!         if (n2p(k,i) < 0.01) n2p(k,i) = 0.01
!
! Convergence Criteria for stopping interation loop (Ne-only)
!          
	  if (it > 1) then
	    test = ABS(1.0-(neold/ne(k,i)))*100.
            ne(k,i)=(ne(k,i)+neold)/2.0
	    if (test < 0.5) go to 100
          endif
!
! Store previous iterated Ne value for next convergence test
!
	  neold = ne(k,i)
!
          enddo ! it=1,100
  100 continue
! --------------------------------------------------------------------

        enddo ! k=lev0,lev1-1
!
        n2p_prod (lev1,i)=spval  ; n2p_loss (lev1,i)=spval
        nop_prod (lev1,i)=spval  ; nop_loss (lev1,i)=spval
        co2p_prod(lev1,i)=spval  ; co2p_loss(lev1,i)=spval
        o2p_prod (lev1,i)=spval  ; o2p_loss (lev1,i)=spval
        op_prod  (lev1,i)=spval  ; op_loss  (lev1,i)=spval
        n2p (lev1,i) = spval
        nop (lev1,i) = spval
        co2p(lev1,i) = spval
        o2p (lev1,i) = spval
        op  (lev1,i) = spval
      enddo ! i=lon0,lon1
!
! -------------------------------------------------------------------
! Convert ion densities to log10 densities for plotting
! Ne and ions (last) and check for values less than 1.0
      do i=lon0,lon1
        do k=lev0,lev1-1
! ------------------------------------------------------------------
	  if (nop(k,i) < 1.0) then
              nopl(k,i) = 0.0
          else
              nopl(k,i)   = log10(nop(k,i))
          endif
! ------------------------------------------------------------------
	  if (o2p(k,i) < 1.0) then
               o2pl(k,i) = 0.0
          else
              o2pl(k,i)   = log10(o2p(k,i))
          endif
! ------------------------------------------------------------------
	  if (op(k,i) < 1.0) then
              opl(k,i) = 0.0
          else
              opl(k,i)   = log10(op(k,i))
          endif
! ------------------------------------------------------------------
	  if (co2p(k,i) < 1.0) then
              co2pl(k,i) = 0.0
          else
              co2pl(k,i)   = log10(co2p(k,i))
          endif
! ------------------------------------------------------------------
	  if (n2p(k,i) < 1.0) then
               n2pl(k,i) = 0.0
          else
               n2pl(k,i)   = log10(n2p(k,i))
          endif
! ------------------------------------------------------------------
	  if (ne(k,i) < 1.0) then
              nel(k,i) = 0.0
          else
              nel(k,i)   = log10(ne(k,i))
          endif
! ------------------------------------------------------------------
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
! -------------------------------------------------------------------
!     call addfld('N2P_P' ,' ',' ',n2p_prod,
!    |  'lev',k0,k1,'lon',i0,i1,lat)
!      call addfld('N2P_L' ,' ',' ',n2p_loss ,
!    |  'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('XN2P',' ',' ',n2pl,'lev',k0,k1,'lon',i0,i1,lat)
!
!     call addfld('NOP_P' ,' ',' ',nop_prod ,
!    |  'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('NOP_L' ,' ',' ',nop_loss ,
!    |  'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('XNOP',' ',' ',nopl,'lev',k0,k1,'lon',i0,i1,lat)
!
!     call addfld('CO2P_P',' ',' ',co2p_prod,
!    |  'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('CO2P_L',' ',' ',co2p_loss,
!    |  'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('XCO2P',' ',' ',co2pl,'lev',k0,k1,'lon',i0,i1,lat)
!
!     call addfld('O2P_P' ,' ',' ',o2p_prod ,
!    |  'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('O2P_L' ,' ',' ',o2p_loss ,
!    |  'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('XO2P',' ',' ',o2pl,'lev',k0,k1,'lon',i0,i1,lat)
!
!     call addfld('OP_P'  ,' ',' ',op_prod  ,
!    |  'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('OP_L'  ,' ',' ',op_loss  ,
!    |  'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('XOP',' ',' ',opl,'lev',k0,k1,'lon',i0,i1,lat)

      call addfld('XNE',' ',' ',nel,'lev',k0,k1,'lon',i0,i1,lat)
! -------------------------------------------------------------------

      end subroutine ionden
      end module ionden_module
!
