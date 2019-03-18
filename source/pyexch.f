      subroutine pyexch(rn,rn13,rnp,rnpp,r,ex,vx)
c-----------------------------------------------------------------------
c     Perdew, Yue          Phys. Rev. B33,8800(1986)
c                          Phys. Rev. B33,8822(1986)
c                          Phys. Rev. B34,7406(1986)
c                 See also Phys. Rev. B43,1399(1991)
c          exchange - contributtion in H
c     Modified by H. Akai, Dec. 2007
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (rm=1d0/15d0,a=0.0864d0/rm,b=14d0,c=2d-1)
      pi=4d0*atan(1d0)
      tpisq13=(3d0*pi**2)**(1d0/3d0)
      ax=-(3d0/4d0)*(3d0/pi)**(1d0/3d0)
      anpp=-rnpp
      ckf=tpisq13*rn13
      s=abs(rnp/(2d0*ckf*rn))
      f0=1d0+a*s**2+b*s**4+c*s**6
      df0=2d0*a*s+4d0*b*s**3+6d0*c*s**5
      f=f0**rm
      ex=ax*rn*rn13*f
      t=(2d0*rnp/r+rnpp)/(2d0*ckf)**2/rn
      u=rnp*anpp/(2d0*ckf)**3/rn**2
      vx=ax*rn13*((4d0/3d0)*f-(t/s)*rm*(f/f0)*df0
     &    -(u-(4d0/3d0)*s**3)*rm*(f/f0)*((rm-1d0)*df0**2/f0/s
     &    +8d0*b*s+24d0*c*s**3))
      end
