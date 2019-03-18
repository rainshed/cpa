      subroutine excev(ro,v,w,dr,xr,msr,mna,meshr)
c-----------------------------------------------------------------------
c     This program generates the generalized gradient approximation for
c     the exchange energy and potential (Engel and Vosko 1993).
c     Energies in Hartrees, distances in Bohrs -> Converted to Ry
c     See: Engel and Vosko, Phys. Rev. B47, 13164 (1993)
c     Routines written by battocletti 19/06/94
c     Modified by H. Akai, Dec. 2007
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 ro(msr,mna,2),v(msr,mna,2),w(msr,mna,*)
     &      ,dr(msr,mna),xr(msr,mna)
      thrd=1d0/3d0
      thrd2=2d0/3d0
      pi=4d0*atan(1d0)
      conf=(3d0*pi**2)**thrd
      conrs=(3d0/4d0/pi)**thrd
      fourpi=4d0*pi
      do 10 ia=1,mna
      do 10 k=1,msr
      exc=0d0
      vxcup=0d0
      vxcdn=0d0
      if(ro(k,ia,1) .gt. 1d-18) then
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
      do 20 isp=1,2
      if(isp .eq. 1 ) then
      d=2d0*rhou
      dp=2d0*rhou1
      dpp=2d0*rhou2
      else
      d=2d0*rhod
      dp=2d0*rhod1
      dpp=2d0*rhod2
      endif
      if(k .eq. meshr)then
      dp=0d0
      dpp=0d0
      endif
      fk=conf*d**thrd
      ss=abs(dp)/(d*2d0*fk)
c     --- xi is the derivative of xi, where xi is given by eq.(4)
c         of Engel and Vosko's.
      dxi=5d-1*dp*d**(-8d0/3d0)*(dpp-(8d0/3d0)*dp**2/d)/conf**2
      uu=dp*dxi/4d0/d/fk**2
      vv=(2d0*dp/xr(k,ia)+dpp)/(d*(2d0*fk)**2)
      call exgga(d,ss,uu,vv,fk,ex,vx)
      exc=exc+ex*(d/2d0)/rho
      if(isp .eq. 1) vxcup=vx
   20 if(isp .eq. 2) vxcdn=vx
c     --- local correlation
      d=rho
      zet=(rhou-rhod)/rho
      rs=conrs/d**thrd
      call corlsd(rs,zet,ec,vcup,vcdn,ecrs,eczet,alfc)
      pi=4d0*atan(1d0)
      exc=exc+ec
      vxcup=vxcup+vcup
      vxcdn=vxcdn+vcdn
c --- nonlocal correlation
      fk=1.91915829d0/rs
      sk=sqrt(4d0*fk/pi)
      g=((1d0+zet)**thrd2+(1d0-zet)**thrd2)/2d0
      tt=abs(rho1)/(d*2d0*sk*g)
      uu=rho1*(-rho2)/(d**2*(2d0*sk*g)**3)
      vv=(2d0*rho1/xr(k,ia)+rho2)/(d*(2d0*sk*g)**2)
      ww=rho1*(rhou1-rhod1-zet*rho1)/(d**2*(2d0*sk*g)**2)
      call corgga(rs,zet,tt,uu,vv,ww,h,dvcup,dvcdn,fk,sk,g,ec
     &           ,ecrs,eczet)
      exc=exc+h
      vxcup=vxcup+dvcup
      vxcdn=vxcdn+dvcdn
      endif
c --- convert from H to Ry
      exc=2d0*exc
      xu=2d0*vxcup
      xd=2d0*vxcdn
      v(k,ia,1)=v(k,ia,1)+xu
      v(k,ia,2)=v(k,ia,2)+xd
      u=w(k,ia,1)
      w(k,ia,1)=ro(k,ia,1)*(u+exc)
c     w(k,ia,2)=ro(k,ia,1)*u-3d0*(exc-xu)*(ro(k,ia,1)+
c    &   ro(k,ia,2))*5d-1-3d0*(exc-xd)*(ro(k,ia,1)-ro(k,ia,2))*5d-1
   10 continue
      end
