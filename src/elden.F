!
      subroutine elden(tn,barm,o1,co,co2,n4s,n2,
     |  n2p,nop,co2p,o2p,oplus,electrons,lev0,lev1,lon0,lon1,lat)
!
! 11/3/06:
! Solve quartic for electron density, and calculate 4 ion species 
! for venus: n2+, no+, co2+, and o2+ (o+ will be added later)
!
      use params_module,only: spval
      use qrj_module,only: qn2p,qco2p,qo2p, ! qnop not available 
     |  mk_xnmbari
      use fields_module,only: tlbc
      use cons_module,only: rmassinv_o1,rmassinv_co,rmassinv_co2,
     |  rmassinv_n4s,rmassinv_n2
      use addfld_module,only: addfld
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
     |  oplus,     ! O+   ! not available yet as of 11/3/06 
     |  electrons  ! electron density (output to f4d(ne))
!
! Local:
      integer :: k,i,k0,k1,i0,i1
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  xo,xco,xco2,xn4s,xn2,  ! number densities
     |  xnmbari,               ! p0*e(-z)*barm/kT at interfaces
     |  tni                    ! tn at interfaces
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  a0,a1,a2,a3,a4, ! coefficients for quartic solver
     |  a,b,c,d,e,fg,h, ! terms for quartic coefficients
     |  root            ! output from quartic solver
!
! For addfld calls:
      k0=lev0 ; k1=lev1 ; i0=lon0 ; i1=lon1
!
! tn at interfaces:
      do i=lon0,lon1
        tni(lev0,i) = tlbc(i,lat)
        do k=lev0+1,lev1-1
          tni(k,i) = .5*(tn(k-1,i)+tn(k,i))
        enddo
        tni(lev1,i) = tn(lev1-1,i) ! nlevp1 <- nlev
      enddo
!     write(6,"('elden: lat=',
!
! Calculate p0*e(-z)*barm/kT (1/(N*MBAR))at interfaces
!     subroutine mk_xnmbari(tni,barm,xnmbari,lev0,lev1,lon0,lon1,lat)
!
      call mk_xnmbari(tni,barm(:,lon0:lon1),xnmbari,
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
!     call addfld('XO'  ,' ',' ',xo  ,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XCO' ,' ',' ',xco ,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XCO2',' ',' ',xco2,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XN4S',' ',' ',xn4s,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XN2' ,' ',' ',xn2 ,'lev',k0,k1,'lon',i0,i1,lat)
!
! Set up terms for quartic coefficients:
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          a(k,i) = 0.

          b(k,i) = 0.

          c(k,i) = 0.

          d(k,i) = 0.

          e(k,i) = 0.

          fg(k,i) = 0.

          h(k,i) = 0.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Coefficients for quartic solver: a0,a1,a2,a3,a4
!
          a0(k,i) = 0.

          a1(k,i) = 0.

          a2(k,i) = 0.

          a3(k,i) = 0.

          a4(k,i) = 0.

        enddo ! k=lev0,lev1-1
        a(lev1,i) = spval
        b(lev1,i) = spval
        c(lev1,i) = spval
        d(lev1,i) = spval
        e(lev1,i) = spval
        fg(lev1,i) = spval
        h(lev1,i) = spval
        a0(lev1,i) = spval
        a1(lev1,i) = spval
        a2(lev1,i) = spval
        a3(lev1,i) = spval
        a4(lev1,i) = spval
      enddo ! i=lon0,lon1
      call addfld('A_COEF' ,' ',' ',a  ,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('B_COEF' ,' ',' ',b  ,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('C_COEF' ,' ',' ',c  ,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('D_COEF' ,' ',' ',d  ,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('E_COEF' ,' ',' ',e  ,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('FG_COEF',' ',' ',fg ,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('H_COEF' ,' ',' ',h  ,'lev',k0,k1,'lon',i0,i1,lat)
!
      call addfld('A0',' ',' ',a0,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('A1',' ',' ',a1,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('A2',' ',' ',a2,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('A3',' ',' ',a3,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('A4',' ',' ',a4,'lev',k0,k1,'lon',i0,i1,lat)
!
! Solve quartic. Vquart returns electron density Ne in root:
!
      call vquart(a0,a1,a2,a3,a4,root,lev0,lev1,lon0,lon1,lat)
!
! Insure positive Ne (at least 3100):
      do i=lon0,lon1
        do k=lev0,lev1-1
          if (root(k,i) < 3.1e3) root(k,i) = 3.1e3
        enddo
      enddo
!
! Calculate ions n2p, nop, co2p, o2p, oplus:
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          n2p  (k,i) = 0.

          nop  (k,i) = 0.

          co2p (k,i) = 0.

          o2p  (k,i) = 0.

          oplus(k,i) = 0.

        enddo ! k=lev0,lev1-1
        n2p  (lev1,i) = spval
        nop  (lev1,i) = spval
        co2p (lev1,i) = spval
        o2p  (lev1,i) = spval
        oplus(lev1,i) = spval
      enddo ! i=lon0,lon1
      call addfld('N2P'  ,' ',' ',n2p  ,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('NOP'  ,' ',' ',nop  ,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('CO2P' ,' ',' ',co2p ,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('O2P'  ,' ',' ',o2p  ,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('OPLUS',' ',' ',oplus,'lev',k0,k1,'lon',i0,i1,lat)
!
! Transfer root to electrons output array:
      do i=lon0,lon1
        do k=lev0,lev1-2
          electrons(k+1,i) = sqrt(root(k,i)*root(k+1,i))
        enddo ! k=lev0,lev1-2
!
! Lower and upper boundaries:
        electrons(lev0,i) = sqrt(root(lev0  ,i)**3/root(lev0+1,i))
        electrons(lev1,i) = sqrt(root(lev1-1,i)**3/root(lev1-2,i))
      enddo ! i=lon0,lon1

      call addfld('NE_ELDEN',' ',' ',electrons(:,i0:i1),
     |  'lev',k0,k1,'lon',i0,i1,lat)

      end subroutine elden
!-----------------------------------------------------------------------
      subroutine vquart(a0,a1,a2,a3,a4,root,lev0,lev1,lon0,lon1,lat)
      implicit none
!
! Determines five roots of the equation:
!   a4*x**4 + 4.*a3*x**3 + 6.*a2*x**2 + 4.*a1*x + a0 = 0.
!
! Procedure is specificlly designed for real quartics with real roots
! only one of which is positive.
!
! This is called by elden for electron density.
! 
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0:lon1),intent(in)  :: a0,a1,a2,a3,a4
      real,dimension(lev0:lev1,lon0:lon1),intent(out) :: root
!
! Local:
      integer :: k,i,nlevs
      real,dimension(lev0:lev1,lon0:lon1) :: w1,w2,w3 ! work arrays
      real,parameter :: e=1.e-300 ! largest exponent on ieee is about 307
!
      nlevs = lev1-lev0+1
      do i=lon0,lon1
        do k=lev0,lev1-1
!
! w1 = ch
          w1(k,i) = -(a4(k,i)*a0(k,i)-4.*a3(k,i)*a1(k,i)+3.*a2(k,i)**2)/
     |                12.
!
! w2 = cg
          w2(k,i) = (a4(k,i)*(a2(k,i)*a0(k,i)-a1(k,i)**2)-a3(k,i)*
     |      (a3(k,i)*a0(k,i)-a1(k,i)*a2(k,i))+a2(k,i)*(a3(k,i)*a1(k,i)-
     |      a2(k,i)**2))/4.
!
! root=rlam=-2.*real((.5*(cmplx(cg,0.)+csqrt(cmplx(cg**2+4.
!      *ch**3+e,0.)))+cmplx(e,0.))**(1./3.))
! 
          root(k,i) = -2.*real((.5*(cmplx(w2(k,i),0.)+
     |      csqrt(cmplx(w2(k,i)**2+4.*w1(k,i)**3+e,0.)))+
     |      cmplx(e,0.))**(1./3.))
!
! W1=P=SQRT(A(5)*RLAM+A(4)**2-A(5)*A(3)+E)
!
          w1(k,i) = a4(k,i)*root(k,i)+a3(k,i)**2-a4(k,i)*a2(k,i)+e
          if (w1(k,i) < 0.) w1(k,i) = 0.
          w1(k,i) = sqrt(w1(k,i))
!
! W2=Q=SQRT((2.*RLAM+A(3))**2-A(5)*A(1)+E)
!
          w2(k,i) = sqrt((2.*root(k,i)+a2(k,i))**2-a4(k,i)*a0(k,i)+e)
!
! W3=PQ=2.*A(4)*RLAM+A(4)*A(3)-A(5)*A(2)+E
!
          w3(k,i) = 2.*a3(k,i)*root(k,i)+a3(k,i)*a2(k,i)-a4(k,i)*a1(k,i)
     |      +e
!
!  W1=P=SIGN(P,Q*PQ)
!
          w1(k,i) = sign(w1(k,i),w2(k,i)*w3(k,i))
!
! W3=P-A4
!
          w3(k,i) = w1(k,i)-a3(k,i)
!
! Final evaluation of root:
!
          root(k,i) = (w3(k,i)+sqrt(w3(k,i)**2-a4(k,i)*(a2(k,i)+2.*
     |      root(k,i)-w2(k,i))))/a4(k,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfsech('W1'    ,' ',' ',w1,lon0,lon1,nlevs,nlevs-1,lat)
!     call addfsech('W2'    ,' ',' ',w2,lon0,lon1,nlevs,nlevs-1,lat)
!     call addfsech('W3'    ,' ',' ',w3,lon0,lon1,nlevs,nlevs-1,lat)
!     call addfsech('VQROOT',' ',' ',root,lon0,lon1,nlevs,nlevs-1,lat)

      end subroutine vquart
