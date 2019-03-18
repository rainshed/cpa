      subroutine guessz(anclr,z,xr,msr,ns1)
c-----------------------------------------------------------------------
c     Gives starting effective nuclear charge for atomic calculations.
c     coded by H.Akai, 1986, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension z(ns1),xr(ns1)
      data pi/3.1415926535898d0/
      do 10 k=2,msr
   10 z(k)=tfp(xr(k),anclr)
      do 20 k=2,msr
      ro=-4d0*z(k)*sqrt(-z(k))/3d0/pi
      z(k)=-(z(k)+fldf(ro))*xr(k)
   20 if(z(k) .lt. 2d0) z(k)=2d0
      z(1)=2d0*anclr
      return
      end
