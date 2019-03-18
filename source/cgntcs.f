      subroutine cgntcs(e,ew,ez,tc,ng)
c-----------------------------------------------------------------------
c     Generate thebycheff polinomials for a complex argument e.
c     coded by H.Akai, 1983, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension tc(ng)
      complex*16 e,tc
      tc(1)=1d0
      tc(2)=(e-ew)/ez
      do 10 n=3,ng
   10 tc(n)=2d0*tc(2)*tc(n-1)-tc(n-2)
      return
      end
