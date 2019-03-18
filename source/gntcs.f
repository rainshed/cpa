      subroutine gntcs(e,ew,ez,tc,ng)
c-----------------------------------------------------------------------
c     Generate Tchebycheff polinomials for a real argument e.
c     coded by M.Akai, 1980, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension tc(ng)
      tc(1)=1d0
      tc(2)=(e-ew)/ez
      do 10 n=3,ng
   10 tc(n)=2d0*tc(2)*tc(n-1)-tc(n-2)
      return
      end
