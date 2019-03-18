      subroutine gcor(a,a1,b1,b2,b3,b4,p,rs,gg,ggrs)
c-----------------------------------------------------------------------
c     called by subroutine corlsd
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      p1=p+1d0
      q0=-2d0*a*(1d0+a1*rs)
      rs12=sqrt(rs)
      rs32=rs12**3
      rsp=rs**p
      q1=2d0*a*(b1*rs12+b2*rs+b3*rs32+b4*rs*rsp)
      q2=dlog(1d0+1d0/q1)
      gg=q0*q2
      q3=a*(b1/rs12+2d0*b2+3d0*b3*rs12+2d0*b4*p1*rsp)
      ggrs=-2d0*a*a1*q2-q0*q3/(q1**2+q1)
      end
