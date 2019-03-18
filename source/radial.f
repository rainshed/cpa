      subroutine radial(e,j,g,p,q,v,w,dr,xr,meshr,isr)
c-----------------------------------------------------------------------
c     Soleves the Schroedinger equation. Suitable for valence
c     states.
c     coded by H.Akai, 1983, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (c1=274.0720442d0, c2=1d0/c1**2, c3=24d0/9d0
     &          ,c4=19d0/9d0, c5=5d0/9d0, c6=1d0/9d0 )
      dimension v(meshr),dr(meshr),xr(meshr),g(meshr),dp(3),dq(3)
     &          ,w(meshr)
      logical sra
      sra=isr .eq. 1
      l=j-1
      bl=dble(l*j)
      es=v(meshr)+e
c
c     ---- expansion of type p=x**(l+1)+...  ----
      do 610 k=1,3
      ek2=(w(k)-es)*xr(k)
      ek3=(v(k)-es)*xr(k)
      rm=xr(k)
      if(sra) rm=rm-ek3*c2
      ek2=bl/rm+ek2
      delt=dr(k)/xr(k)
      p=xr(k)**j
      q=dble(l)*p/rm
      dp(k)=(p+rm*q)*delt/c3
      dq(k)=(ek2*p-q)*delt/c3
      g(k)=p
c     --- if the small component is to be take into account,
c         the following line must be activatd.
  610 if(sra) g(k)=g(k)*sqrt(1d0+bl/(rm*c1)**2+(q/p/c1)**2)
      s=5d-1*g(1)**2*dr(1)+2d0*g(2)**2*dr(2)+g(3)**2*dr(3)
c
c     ---- adams-multon ----
      do 20 k=4,meshr-1
      ek2=(w(k)-es)*xr(k)
      ek3=(v(k)-es)*xr(k)
      rm=xr(k)
      if(sra) rm=rm-ek3*c2
      ek2=bl/rm+ek2
      dlt=c3*xr(k)/dr(k)
      gm=(1d0-dlt)/ek2
      pm=p+c4*dp(3)-c5*dp(2)+c6*dp(1)
      qm=q+c4*dq(3)-c5*dq(2)+c6*dq(1)
      q=(gm*qm-pm)/(rm+gm*(dlt+1d0))
      p=(pm+rm*q)/(dlt-1d0)
      dp(1)=dp(2)
      dp(2)=dp(3)
      dq(1)=dq(2)
      dq(2)=dq(3)
      dp(3)=p+rm*q
      dq(3)=ek2*p-q
      p=p*dlt
      q=q*dlt
      g(k)=p
c     --- if the small component is to be take into account,
c         the following line must be activatd.
      if(sra) g(k)=g(k)*sqrt(1d0+bl/(rm*c1)**2+(q/p/c1)**2)
      ps=dr(k)*g(k)**2
   20 s=s+dble(mod(k+1,2)+1)*ps
c
c
      s=1d0/sqrt((2d0*s-ps)/3d0)
      do 30 k=1,meshr-1
   30 g(k)=g(k)*s
      g(meshr)=0d0
      q=s*(dp(3)*dlt-p)/xr(meshr-1)**2
      p=p*s/xr(meshr-1)
c     write(6,2000)e,l,g(1),g(2),g(3),g(4),g(5),g(6)
c2000 format(' e=',f5.1,' l=',i1,' r=',1p,6e14.6)
      end
