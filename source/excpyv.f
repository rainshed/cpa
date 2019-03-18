      subroutine excpyv(ro,v,w,dr,xr,msr,mna,meshr)
c-----------------------------------------------------------------------
c     Perdew, Yue          Phys. Rev. B33,8800(1986)
c                          Phys. Rev. B33,8822(1986)
c                          Phys. Rev. B34,7406(1986)
c                 see also Phys. Rev. B43,1399(1991)
c          sdf parametrization by vwn
c     Modified by H. Akai, Dec. 2007
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 ro(msr,mna,2),v(msr,mna,2),w(msr,mna,*)
     &      ,dr(msr,mna),xr(msr,mna)
      real*8 zeta, n13, nu13,nd13
      data cg1/0.001667d0/, cg2/0.002568d0/, cga/0.023266d0/
     &    ,cgb/7.389d-6/, cgg/8.723d0/, cgd/0.472d0/
      t=1d0/3d0
      pi=4d0*atan(1d0)
      ftl=1.1d-1
      bpy=1.745415107d0*ftl
      vfd=sqrt(5d-1)
      t23=2d0**(2d0/3d0)
      t13=2d0**t
      p73=7d0/3d0
      p76=7d0/6d0
      cninf=cg1+cg2
      fourpi=4d0*pi
      do 10 ia=1,mna
      do 10 k=1,msr
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
      call vwncor(rs,rt,xu,xd,exc,vculoc,vcdloc,ecloc)
c     --- non-local corrections
      if(k .le. meshr-1) then
      n13=rho**t
      nu13=rhou**t
      nd13=rhod**t
      zeta=(rhou-rhod)/rho
      if(zeta .gt. 1d0) zeta=1d0
      if(zeta .lt. -1d0) zeta=-1d0
      d=vfd*sqrt((1d0+zeta)**(5d0/3d0)
     & +(1d0-zeta)**(5d0/3d0))
      call pyexch(2d0*rhou,t13*nu13,2d0*rhou1,2*rhou2,xr(k,ia),exu,vxu)
      call pyexch(2d0*rhod,t13*nd13,2d0*rhod1,2*rhod2,xr(k,ia),exd,vxd)
      exnl=5d-1*(exu+exd)/rho
      cnz=cg2+cga*rs+cgb*rs**2
      cnn=1d0+cgg*rs+cgd*rs**2+1d4*cgb*rs**3
      cn=cg1+cnz/cnn
      fpy=bpy*(cninf/cn)*abs(rho1)/rho**p76
      epy=exp(-fpy)
      ecnl=(1d0/d)*epy*cn*abs(rho1)**2/rho**p73
      excnl=exnl+ecnl
      xn=rho1/rho
      yn=(2d0*rho1/xr(k,ia)+rho2) /rho
      danr=-rho2
      p2=(2d0-fpy)*yn
      p3=(4d0/3d0-11d0*fpy/3d0+7d0*fpy*fpy/6d0)*xn**2
      p4=fpy*(fpy-3d0)*rho1*danr/rho/abs(rho1)
      p5u=5d0*n13*(rhou**(2d0/3d0)-rhod**(2d0/3d0))
     &   *rho1/(6d0*d**2*rho**4)
      p5d=-p5u
      cnzp=cga+2d0*cgb*rs
      cnnp=cgg+2d0*cgd*rs+3d0*1d4*cgb*rs**2
      p7=(rho1**2/rho)*(fpy**2-fpy-1d0)
     &  *(cnzp-cnz*cnnp/cnn)/cnn*(-4d0*pi/9d0)*rs**4
      p6u=(1d0-fpy)*rhod*rho1-(2d0-fpy)*rho*rhod1
      vnlu=vxu-epy/d/n13*(cn*(p2-p3+p4-p5u*t23*p6u)-p7)
      p6d=(1d0-fpy)*rhou*rho1-(2d0-fpy)*rho*rhou1
      vnld=vxd-epy/d/n13*(cn*(p2-p3+p4-p5d*t23*p6d)-p7)
c     --- convert nl-contributions from H to Ry
      exc=ecloc+2d0*excnl
      xu=vculoc+2d0*vnlu
      xd=vcdloc+2d0*vnld
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
c     w(k,ia,2)=ro(k,ia,1)*u-3d0*(exc-xu)*(ro(k,ia,1)+ro(k,ia,2))*5d-1
c    &                      -3d0*(exc-xd)*(ro(k,ia,1)-ro(k,ia,2))*5d-1
   10 continue
      end
