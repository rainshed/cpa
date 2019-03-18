      function dawson(x)
c-----------------------------------------------------------------------
c     given x, this function returns Dawson's integral F(x) defined by
c       F(x)=exp(-x^2)*integral(0 through x) exp(t^2) dt
c     with fractional error less than 1d-14.
c     
c     coded by H. Akai 9 April 2002, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(nmx1=12,nmx2=16, h=2d-1)
      real*8 a(nmx1),b(nmx2)
      logical init
      save init,a,b,ci
      data init/.true./
      if(init) then
      init=.false.
      pi=4d0*atan(1d0)
      ci=1d0/sqrt(pi)
      do 10 i=1,nmx1
   10 a(i)=-2d0/dble(i+i+1)
      do 20 i=1,nmx2
   20 b(i)=exp(-(dble(i+i-1)*h)**2)
      endif
      if(abs(x) .lt. 6d-1) then
c---  series expansion is accurate for abs(x)<0.6
      xx=x**2
      dawson=1d0
      do 30 i=nmx1,1,-1
   30 dawson=dawson*a(i)*xx+1d0
      dawson=dawson*x
      return
      else
c---  Rybicki's method (see Numerical Receipes) is used.
      n0=0.5d0*abs(x)/h
      n0=n0+n0
      s=abs(x)-dble(n0)*h
      cen=exp(2d0*s*h)
      ce=cen**2
      fct=sign(exp(-s**2),x)
      dawson=0d0
      do 40 i=1,nmx2
      nn=i+i-1
      dawson=dawson+b(i)*(cen/dble(n0+nn)+1d0/dble(n0-nn)/cen)
   40 cen=cen*ce
      dawson=dawson*ci*fct
      endif
      end
