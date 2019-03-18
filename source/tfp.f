      function tfp(r,z)
c-----------------------------------------------------------------------
c     Thomas-Fermi potential   w=z-numbre of electrons
c     coded by H.Akai, 1985, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      data w/0d0/
      x=sqrt((r*(z+w)**(1d0/3d0))/0.8853d0)
      d=x*(0.60112d0*x+1.81061d0)+1d0
      e=x*(x*(x*(x*(x*0.04793d0+0.21465d0)+0.77112d0)
     &   +1.39515d0)+1.81061d0)+1d0
      ze=(z+w)*(d/e)**2-w
      tfp=-2d0*ze/r
      return
      end
