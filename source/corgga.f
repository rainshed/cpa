      subroutine corgga(rs,zet,t,uu,vv,ww,h,dvcup,dvcdn,fk,sk,g,ec
     &                 ,ecrs,eczet)
c-----------------------------------------------------------------------
c     gga91 correlation
c     input rs: seitz radius
c     input zet: relative spin polarization
c     input t: abs(grad d)/(d*2.*ks*g)
c     input uu: (grad d)*grad(abs(grad d))/(d**2 * (2*ks*g)**3)
c     input vv: (laplacian d)/(d * (2*ks*g)**2)
c     input ww:  (grad d)*(grad zet)/(d * (2*ks*g)**2
c     output h: nonlocal part of correlation energy per electron
c     output dvcup,dvcdn:  nonlocal parts of correlation potentials
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      data xnu,cc0,cx,alf/15.75592d0,0.004235d0,-0.001667212d0,0.09d0/
     &    ,c1,c2,c3,c4/0.002568d0,0.023266d0,7.389d-6,8.723d0/
     &    ,c5,c6,a4/0.472d0,7.389d-2,1d2/, zero/1d-10/
      if(abs(abs(zet)-1d0) .lt. zero) zet=sign(1d0-zero,zet)
      thrdm=-1d0/3d0
      thrd2=2d0/3d0
      bet=xnu*cc0
      delt=2d0*alf/bet
      g3=g**3
      g4=g3*g
      pon=-delt*ec/(g3*bet)
      b=delt/(dexp(pon)-1d0)
      b2=b**2
      t2=t**2
      t4=t2**2
      t6=t2**3
      rs2=rs*rs
      rs3=rs2*rs
      q4=1d0+b*t2
      q5=1d0+b*t2+b2*t4
      q6=c1+c2*rs+c3*rs2
      q7=1d0+c4*rs+c5*rs2+c6*rs3
      cc=-cx+q6/q7
      r0=(sk/fk)**2
      r1=a4*r0*g4
      coeff=cc-cc0-3d0*cx/7d0
      r2=xnu*coeff*g3
      r3=dexp(-r1*t2)
      h0=g3*(bet/delt)*log(1d0+delt*q4*t2/q5)
      h1=r3*r2*t2
      h=h0+h1
c     ---local correlation option:
c      h=0d0
c     ---energy done. now the potential:
      ccrs=(c2+2d0*c3*rs)/q7-q6*(c4+2d0*c5*rs+3d0*c6*rs2)/q7**2
      rsthrd=rs/3d0
      r4=rsthrd*ccrs/coeff
      gz=((1d0+zet)**thrdm-(1d0-zet)**thrdm)/3d0
      fac=delt/b+1d0
      bg=-3d0*b2*ec*fac/(bet*g4)
      bec=b2*fac/(bet*g3)
      q8=q5*q5+delt*q4*q5*t2
      q9=1d0+2d0*b*t2
      h0b=-bet*g3*b*t6*(2d0+b*t2)/q8
      h0rs=-rsthrd*h0b*bec*ecrs
      fact0=2d0*delt-6d0*b
      fact1=q5*q9+q4*q9*q9
      h0bt=2d0*bet*g3*t4*((q4*q5*fact0-delt*fact1)/q8)/q8
      h0rst=rsthrd*t2*h0bt*bec*ecrs
      h0z=3d0*gz*h0/g+h0b*(bg*gz+bec*eczet)
      h0t=2d0*bet*g3*q9/q8
      h0zt=3d0*gz*h0t/g+h0bt*(bg*gz+bec*eczet)
      fact2=q4*q5+b*t2*(q4*q9+q5)
      fact3=2d0*b*q5*q9+delt*fact2
      h0tt=4d0*bet*g3*t*(2d0*b/q8-(q9*fact3/q8)/q8)
      h1rs=r3*r2*t2*(-r4+r1*t2/3d0)
      fact4=2d0-r1*t2
      h1rst=r3*r2*t2*(2d0*r4*(1d0-r1*t2)-thrd2*r1*t2*fact4)
      h1z=gz*r3*r2*t2*(3d0-4d0*r1*t2)/g
      h1t=2d0*r3*r2*(1.d0-r1*t2)
      h1zt=2d0*gz*r3*r2*(3d0-11d0*r1*t2+4d0*r1*r1*t4)/g
      h1tt=4d0*r3*r2*r1*t*(-2d0+r1*t2)
      hrs=h0rs+h1rs
      hrst=h0rst+h1rst
      ht=h0t+h1t
      htt=h0tt+h1tt
      hz=h0z+h1z
      hzt=h0zt+h1zt
      comm=h+hrs+hrst+t2*ht/6d0+7d0*t2*t*htt/6d0
      pref=hz-gz*t2*ht/g
      fact5=gz*(2d0*ht+t*htt)/g
      comm=comm-pref*zet-uu*htt-vv*ht-ww*(hzt-fact5)
      dvcup=comm+pref
      dvcdn=comm-pref
c     ---local correlation option:
c      dvcup=0d0
c      dvcdn=0d0
c
      end
