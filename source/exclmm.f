      subroutine exclmm(ro,v,w,dr,xr,msr,mna,meshr)
c----------------------------------------------------------------------
c     Langreth-Perdew-Mehl Phys. Rev. B28,1809(1983)
c     See Kutzler Painter  Phys. Rev. B37,2850(1988)
c     sdf parametrization by mjw
c     Adapted to kkr-code
c     Modified by H. Akai, Dec. 2007
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 ro(msr,mna,2),v(msr,mna,2),w(msr,mna,*)
     &      ,dr(msr,mna),xr(msr,mna)
      data a/7.93700526d-1/, gamma/5.129762803d0/
     &    ,cp,rp,cf,rf/4.5d-2,  21d0,  2.25d-2,  52.916684096d0/
c    &    ,cp,rp,cf,rf/5.04d-2, 30d0,  2.54d-2,     75d0/
      g(x)=log(1d0+1d0/x)
      f(c,x)=((c*x-1d0)*x+5d-1)*x+c-t
      pi=4d0*atan(1d0)
      t=1d0/3d0
c     ---     a=2d0**(-t)
c             gamma=(4d0/3d0)*a/(1d0-a)
c         (  for mjw parameter  cf=cp/2 and rf=2**(4/3)*a/(1-a) holds )
c
c     --- alm: 2* the value of  kutzler painter to get ry instead of h
      fcut=1.5d-1
      alm=4.287115843d-03
      blm=1.745415107d0*fcut
      vfd=sqrt(5d-1)
      vf1=-7d0/(9d0*2d0**t)
      vf6=2d0**(2d0/3d0)
      fourpi=4d0*pi
      do 10 ia=1,mna
      do 10 k=1,meshr
      if(ro(k,ia,1) .gt. 1d-20 .and.
     &       .not. abs(ro(k,ia,2)) .gt. ro(k,ia,1)) then
c     --- preparation
      rho=ro(k,ia,1)/fourpi
      sho=ro(k,ia,2)/fourpi
      if(k .le. meshr-1) then
      call drivtv(ro(1,ia,1),meshr-1,xr(1,ia),dr(1,ia),k,rho1,rho2)
      call drivtv(ro(1,ia,2),meshr-1,xr(1,ia),dr(1,ia),k,sho1,sho2)
      else
      rho1=0d0
      rho2=0d0
      sho1=0d0
      sho2=0d0
      endif
      rho1=rho1/fourpi
      rho2=rho2/fourpi
      sho1=sho1/fourpi
      sho2=sho2/fourpi
      rhou=5d-1*(rho+sho)
      rhod=5d-1*(rho-sho)
      rhou1=5d-1*(rho1+sho1)
      rhod1=5d-1*(rho1-sho1)
      rhou2=5d-1*(rho2+sho2)
      rhod2=5d-1*(rho2-sho2)
c     --- main part starts here
      rs=(3d0/ro(k,ia,1))**t
      rt=ro(k,ia,2)/ro(k,ia,1)
      rt=min(1d0,max(-1d0,rt))
      s=5d-1*(1d0+rt)
      ff=(s**(4d0/3d0)+(1d0-s)**(4d0/3d0)-a)/(1d0-a)
      x=rs/rp
      y=rs/rf
      fxp=-1.221774115422d0/rs
      eexp=fxp*.75d0
      cpl=g(x)
      cfl=g(y)
      fcp=-cp*cpl
      fcf=-cf*cfl
      ecp=-cp*f(cpl,x)
      ecf=-cf*f(cfl,y)
      cnu=gamma*(ecf-ecp)
      exc=eexp+ecp+(fxp+cnu)*ff/gamma
      c1=fxp+cnu
      c2=fcp-cnu
      tcf=(fcf-fcp-4d0*(ecf-ecp)/3d0)*ff
      xu=c1*(1d0+rt)**t+c2+tcf
      xd=c1*(1d0-rt)**t+c2+tcf
c     --- non-local corrections
      if(k .le. meshr-1) then
      flm=blm*abs(rho1)/rho**(7d0/6d0)
      elm=2d0*exp(-flm)
      zeta=(rhou-rhod)/rho
      if(zeta .gt. 1d0) zeta=1d0
      if(zeta .lt. -1d0) zeta=-1d0
      d=vfd*sqrt((1d0+zeta)**(5d0/3d0)
     &  +(1d0-zeta)**(5d0/3d0) )
      excnl=alm*(vf1*(rhou1**2/rhou**(4d0/3d0)
     &     +rhod1**2/rhod**(4d0/3d0))/rho
     &     +(elm/d)*rho1**2/rho**(7d0/3d0))
      xn=rho1/rho
      yn=(2d0*rho1/xr(k,ia)+rho2)/rho
      rhothrd=rho**t
      danr=-rho2
      p2=(2d0-flm)*yn
      p3=(4d0/3d0-11d0*flm/3d0+7d0*flm*flm/6d0)*xn**2
      p4=flm*(flm-3d0)*rho1*danr/rho/abs(rho1)
      p5u=5d0*rhothrd*(rhou**(2d0/3d0)-rhod**(2d0/3d0))
     &      *rho1/(6d0*d**2*rho**4)
      p5d=-p5u
      xnu=rhou1/rhou
      ynu=(2d0*rhou1/xr(k,ia)+rhou2)/rhou
      p1u=vf1*((4d0/3d0)*xnu**2-2d0*ynu)*(rhou/rho)**(-t)
      p6u=(1d0-flm)*rhod*rho1-(2d0-flm)*rho*rhod1
      vnlu=(alm/rhothrd)*(p1u-(elm/d)*(p2-p3+p4-p5u*vf6*p6u))
      xnd=rhod1/rhod
      ynd=(2d0*rhod1/xr(k,ia)+rhod2)/rhod
      p1d=vf1*((4d0/3d0)*xnd**2-2d0*ynd)*(rhod/rho)**(-t)
      p6d=(1d0-flm)*rhou*rho1-(2d0-flm)*rho*rhou1
      vnld=(alm/rhothrd)*(p1d-(elm/d)*(p2-p3+p4-p5d*vf6*p6d))
      exc=exc+excnl
      xu=xu+vnlu
      xd=xd+vnld
      endif
      else
      exc=0d0
      xu=0d0
      xd=0d0
      endif
      v(k,ia,1)=v(k,ia,1)+xu
      v(k,ia,2)=v(k,ia,2)+xd
      u=w(k,ia,1)
      w(k,ia,1)=ro(k,ia,1)*(u+exc)
c     w(k,ia,2)=ro(k,ia,1)*u
c    &         -3d0*(exc-xu)*(ro(k,ia,1)+ro(k,ia,2))*5d-1
c    &         -3d0*(exc-xd)*(ro(k,ia,1)-ro(k,ia,2))*5d-1
   10 continue
      end
