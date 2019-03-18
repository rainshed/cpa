      subroutine extorg(org,a,x)
c-----------------------------------------------------------------------
c     Extrapolate charge density toward the origin.
c     coded by H.Akai, 1983, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(3),x(3)
      p=((a(1)-a(2))/(x(1)-x(2))-(a(3)-a(2))/(x(3)-x(2)))/(x(1)-x(3))
      q=(a(1)-a(2))/(x(1)-x(2))-p*(x(1)-x(2))
      org=p*x(2)**2-q*x(2)+a(2)
      return
      end
