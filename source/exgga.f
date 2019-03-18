      subroutine exgga(d,s,u,v,fk,ex,vx)
c-----------------------------------------------------------------------
c     gga Eengel & Vosko exchange for a spin-unpolarized electronic
c     system
c     input d : density
c           s:  abs(grad d)/(2*fk*d)
c           u:  grad(abs(grad d))/(d)
c           v: (laplacian d)/(d*(2*fk)**2)
c           fk: fermi momentum
c     output:  exchange energy per electron (ex) and potential (vx)
c     Modified by H. Akai, Dec. 2007
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      data a1,a2,a3/1.647127,0.980118 ,0.017399/
     &    ,b1,b2,b3/1.523671,0.367229,0.011282/
     &    ,ax/-0.7385588d0/
      pi=4d0*atan(1d0)
      thrd=1d0/3d0
      thrd4=4d0/3d0
      fac=ax*d**thrd
      c2=a1-b1
      s2=s**2
      s4=s2**2
      s6=s4*s2
c     --- pade-approximant
      anomin=1d0+a1*s2+a2*s4+a3*s6
      denomi=1d0+b1*s2+b2*s4+b3*s6
      f=anomin/denomi
      ex=fac*f
c     --- local exchange option
      ex=fac
c     --- energy done. now the potential:
c     --- first derivation of f(s2) with respect to s2
      dn=a1+2*a2*s2+3*a3*s4
      dd=b1+2*b2*s2+3*b3*s4
      df=(dn-dd*f)/denomi
c     --- second derivation of f(s2)
      ddn=2d0*a2+6d0*a3*s2
      ddd=2d0*b2+6d0*b3*s2
      ddf=(ddn-(dn*dd/denomi)-ddd*f-dd*df+(dd*dd*f/denomi))/denomi
      vlda=-fk/pi
      tau=u
c     tau=s2*(u-8d0/3d0*s2)/2d0/fk
      vx=vlda*(f-1.5d0*(df*v+ddf*tau))
c     --- local exchange option:
      vx=fac*thrd4
      end
