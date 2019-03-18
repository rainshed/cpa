      subroutine uxcor(rs,s,uxc1,uxc2,exc)
c-----------------------------------------------------------------------
c    ---- This subroutine was coded by M. Mannien ----
c    subroutine ( M. Manninen )  to calculate energy and potential
c    from Ceperley-Alder ( parametrization of Vosko, Wilk and Nusair )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      data ap,xp0,bp,cp,qp,cp1,cp2,cp3/0.0621814d0,-0.10498d0,
     &  3.72744d0,12.9352d0,6.1519908d0,1.2117833d0,1.1435257d0,
     &  -0.031167608d0/
      data af,xf0,bf,cf,qf,cf1,cf2,cf3/0.0310907d0,-0.32500d0,
     &  7.06042d0,18.0578d0,4.7309269d0,2.9847935d0,2.7100059d0,
     &  -0.1446006d0/
      x=sqrt(rs)
      xpx=x*x+bp*x+cp
      xfx=x*x+bf*x+cf
      s4=s**4-1d0
      cbrt1=(1d0+s)**(1d0/3d0)
      cbrt2=(1d0-s)**(1d0/3d0)
      fs=((1d0+s)**(4d0/3d0)+(1d0-s)**(4d0/3d0)-2d0)/
     &  ( 2d0**(4d0/3d0)-2d0)
      beta=1d0/(2.74208d0+3.182d0*x+0.09873d0*x*x+0.18268d0*x**3)
      dfs=4d0/3d0*(cbrt1-cbrt2)/(2d0**(4d0/3d0)-2d0)
      dbeta=-(0.27402d0*x+0.09873+1.591d0/x)*beta**2
      atnp=atan(qp/(2d0*x+bp))
      atnf=atan(qf/(2d0*x+bf))
      ecp=ap*(log(x*x/xpx)+cp1*atnp
     &   -cp3*(log((x-xp0)**2/xpx)+cp2*atnp))
      ecf=af*(log(x*x/xfx)+cf1*atnf
     &   -cf3*(log((x-xf0)**2/xfx)+cf2*atnf))
      ec=ecp+fs*(ecf-ecp)*(1d0+s4*beta)
      tp1=(x*x+bp*x)/xpx
      tf1=(x*x+bf*x)/xfx
      ucp=ecp-ap/3d0*(1d0-tp1-cp3*(x/(x-xp0)-tp1-xp0*x/xpx))
      ucf=ecf-af/3d0*(1d0-tf1-cf3*(x/(x-xf0)-tf1-xf0*x/xfx))
      uc0=ucp+(ucf-ucp)*fs
      uc10=uc0-(ecf-ecp)*(s-1d0)*dfs
      uc20=uc0-(ecf-ecp)*(s+1d0)*dfs
      duc=(ucf-ucp)*beta*s4*fs
     &  +(ecf-ecp)*(-rs/3d0)*dbeta*s4*fs
      duc1=duc-(ecf-ecp)*beta*(s-1d0)*(4d0*s**3*fs+s4*dfs)
      duc2=duc-(ecf-ecp)*beta*(s+1d0)*(4d0*s**3*fs+s4*dfs)
      uc1=uc10+duc1
      uc2=uc20+duc2
      uxc1=uc1-1.221774d0/rs*cbrt1
      uxc2=uc2-1.221774d0/rs*cbrt2
      exc=ec-0.9163306/rs-0.2381735d0/rs*fs
      return
      end
