      function fintgr(f,dr,xr,meshr)
c-----------------------------------------------------------------------
c     Return radial integral of f.
c     coded by H.Akai
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 f(meshr),dr(meshr),xr(meshr)
      s1=(f(1)*dr(1)*xr(1)**2
     &       +f(meshr-1)*dr(meshr-1)*xr(meshr-1)**2)/2d0
      do 10 k=3,meshr-3,2
   10 s1=s1+f(k)*dr(k)*xr(k)**2
      s2=0d0
      do 20 k=2,meshr-2,2
   20 s2=s2+f(k)*dr(k)*xr(k)**2
      fintgr=(2d0*s1+4d0*s2)/3d0
      return
      end
