      subroutine poisnb(ro,z,anclr,a,b,xr,msr,ns1)
c-----------------------------------------------------------------------
c      This program gives the solution of the Poisson equation for a
c      given charge density.  The Noumerov method is used.
c
c       ro      ... 4*pi times charge density
c       z       ... solution v has the form v = -z(x)/x
c       anclr   ... point charge at the origin neutralizing the system
c       a,b,xr  ... radial mesh must have a form xr(s)=a*(exp(b*s)-1)
c       msr     ... s = 1,2,3,...,msr
c       ns1     ... array size for r0,z,xr
c
c      A change of the variables x ---> s and z ---> y=z/sqrt(x')
c      is employed.  With these new variables the original equation
c      becomes y''(s) = b**2*y(s)/4 - 2*ro(s)*x(s)*x'(s)**(3/2)
c      which can be solved by the Noumerov method.
c
c            coded by H.Akai   Jan. 9, 1986   (Osaka)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ro(ns1),z(ns1),xr(ns1)
      z(1)=0d0
      z1=0d0
      z(2)=-xr(2)**2*ro(1)/3d0
      dr=b*(a+xr(2))
      sdr=sqrt(dr)
      z2=z(2)/sdr
      r1=0d0
      r2=ro(2)*xr(2)*dr*sdr
      bb4=b**2/4d0
      c1=2d0/(12d0-bb4)
      c2=c1*(12d0+5d0*bb4)
      do 10 k=3,msr
      dr=b*(a+xr(k))
      sdr=sqrt(dr)
      r3=ro(k)*xr(k)*dr*sdr
      z3=c2*z2-z1-c1*(r3+10d0*r2+r1)
      z(k)=z3*sdr
      z1=z2
      z2=z3
      r1=r2
   10 r2=r3
      z0=2d0*anclr
      c3=(z0-z(msr))/xr(msr)
      do 20 k=1,msr
   20 z(k)=z0-z(k)-c3*xr(k)
      return
      end
