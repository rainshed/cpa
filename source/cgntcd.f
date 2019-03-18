      subroutine cgntcd(e,ew,ez,tc,td,ng)
c-----------------------------------------------------------------------
c     Generate Tchebycheff polinomials and their derivatives for
c     a complex argument e.
c     coded by H.Akai, 1983, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension tc(ng),td(ng)
      complex*16 e,tc,td
      tc(1)=1d0
      tc(2)=(e-ew)/ez
      td(1)=0d0
      td(2)=1d0/ez
      do 10 n=3,ng
      tc(n)=2d0*tc(2)*tc(n-1)-tc(n-2)
   10 td(n)=2d0*tc(2)*td(n-1)+2d0*tc(n-1)/ez-td(n-2)
      return
      end
