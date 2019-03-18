      subroutine vrotat(v2,v1,alpha,beta,gamma)
c-----------------------------------------------------------------------
c     Rotates a vector v1 arond fixed axis by (alpha,beta,gamma),
c     first around z-axis by alpha, then around old y-axis by beta,
c     and finally again around old z-axis by gamma, giving new vector
c     v2. This corresponds to the rotation of the vector by an Euler
c     angle (gamma,beta,alpha), i.e., first around z-axis by gamma,
c     then around new y-axis by beta, and finally around new z-axis
c     by alpha.
c     Coded by H.Akai, July 1992, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 v1(3),v2(3)
c--     (v1(1),v1(2),v1(3)) -> (t,y,v1(3))
      c=cos(alpha)
      s=sin(alpha)
      t=c*v1(1)-s*v1(2)
      y=s*v1(1)+c*v1(2)
c---    (t,y,v1(3)) -> (x,y,v2(3))
      c=cos(beta)
      s=sin(beta)
      v2(3)=c*v1(3)-s*t
      x=s*v1(3)+c*t
c---     (x,y,v2(3)) -> (v2(1),v2(2),v2(3))
      c=cos(gamma)
      s=sin(gamma)
      v2(1)=c*x-s*y
      v2(2)=s*x+c*y
      end
