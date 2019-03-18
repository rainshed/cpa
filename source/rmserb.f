      subroutine rmserb(v1,v2,rms,dr,xr,msr)
c-----------------------------------------------------------------------
c     rms error analysis
c     coded by H.Akai, 1986, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension dr(msr),xr(msr),v1(msr),v2(msr)
      atvol=xr(msr)**3/3d0
      rms=0d0
      ints=msr-1
      do 10 k=2,ints,2
   10 rms=rms+(v1(k)-v2(k))**2*dr(k)
      rms=rms*2d0
      do 20 k=3,ints,2
   20 rms=rms+(v1(k)-v2(k))**2*dr(k)
      rms=rms*2d0+(v1(1)-v2(1))**2*dr(1)
     &    +(v1(msr)-v2(msr))**2*dr(msr)
      rms=sqrt(rms/3d0/atvol)
      return
      end
