      subroutine excgga(d,s,u,v,ex,vx)
c-----------------------------------------------------------------------
c     gga91 exchange for a spin-unpolarized electronic system
c     input d : density
c           s:  abs(grad d)/(2*kf*d)
c           u:  (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
c           v: (laplacian d)/(d*(2*kf)**2)
c     output:  exchange energy per electron (ex) and potential (vx)
c     Modified by H. Akai, Dec. 2007
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      data a1,a2,a3,a4/0.19645d0,0.27430d0,0.15084d0,100.d0/
     &    ,ax,a,b1/-0.7385588d0,7.7956d0,0.004d0/
      thrd=1d0/3d0
      thrd4=4d0/3d0
      fac=ax*d**thrd
      s2=s**2
      s3=s2*s
      s4=s3*s
      p0=1d0/dsqrt(1d0+a**2*s2)
      p1=log(a*s+1d0/p0)
      p2=dexp(-a4*s2)
      p3=1d0/(1d0+a1*s*p1+b1*s4)
      p4=1d0+a1*s*p1+(a2-a3*p2)*s2
      f=p3*p4
      ex=fac*f
c     --- local exchange option
c      ex = fac
c     --- energy done. now the potential:
      p5=b1*s2-(a2-a3*p2)
      p6=a1*s*(p1+a*s*p0)
      p7=2d0*(a2-a3*p2)+2d0*a3*a4*s2*p2-4d0*b1*s2*f
      fs=p3*(p3*p5*p6+p7)
      p8=2d0*s*(b1-a3*a4*p2)
      p9=a1*p1+a*a1*s*p0*(3d0-a*a*s2*p0*p0)
      p10=4d0*a3*a4*s*p2*(2d0-a4*s2)-8d0*b1*s*f-4.d0*b1*s3*fs
      p11=-p3*p3*(a1*p1+a*a1*s*p0+4d0*b1*s3)
      fss=p3*p3*(p5*p9+p6*p8)+2d0*p3*p5*p6*p11+p3*p10+p7*p11
      vx=fac*(thrd4*f-(u-thrd4*s3)*fss-v*fs)
c     --- local exchange option:
c      vx = fac*thrd4
      end
