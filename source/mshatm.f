      subroutine mshatm(a,b,dr,xr,z,msr,ns1)
c--------------------------------------------------------------------
c     Constructing radial mesh suitable for atomic calculations.
c     coded by H.Akai, 1986, Osaka
c--------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension dr(ns1),xr(ns1)
      a=0.0466d0*0.8853d0*z**(-0.3333333333d0)
      b=0.019d0
c     b=0.002d0
      do 10 k=1,msr
      aebx=a*exp(b*dble(k-1))
      xr(k)=aebx-a
   10 dr(k)=aebx*b
      return
      end
