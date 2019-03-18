      subroutine excvbh(ro,v,w,dr,xr,msr,mna,meshr)
c----------------------------------------------------------------------
c     +----------------------------------+
c        SDF parametrization by von B-H
c     +----------------------------------+
c     coded by H.Akai, 1978, Osaka
c     very minor revision by H.Akai, 23 Nov. 1995, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 ro(msr,mna,2),v(msr,mna,2),w(msr,mna,*)
     &      ,dr(msr,mna),xr(msr,mna)
      data t/.333333333333d0/, a/7.93700526d-1/, gamma/5.129762803d0/
c    &    ,cp,rp,cf,rf/4.5d-2,  21d0,  2.25d-2,  52.916684096d0/
     &    ,cp,rp,cf,rf/5.04d-2, 30d0,  2.54d-2,     75d0/
      g(x)=log(1d0+1d0/x)
      f(c,x)=((c*x-1d0)*x+5d-1)*x+c-t
c
c          a=2d0**(-t)
c          gamma=(4d0/3d0)*a/(1d0-a)
c      (  for mjw parameter  cf=cp/2 and rf=2**(4/3)*a/(1-a) holds )
c
      do 10 ia=1,mna
      do 10 k=1,meshr
c
      if(ro(k,ia,1) .gt. 1d-20 .and.
     &       .not. abs(ro(k,ia,2)) .gt. ro(k,ia,1)) then
      rs=(3d0/ro(k,ia,1))**t
      rt=ro(k,ia,2)/ro(k,ia,1)
      rt=min(1d0,max(-1d0,rt))
      s=5d-1*(1d0+rt)
      ff=(s**(4d0/3d0)+(1d0-s)**(4d0/3d0)-a)/(1d0-a)
      x=rs/rp
      y=rs/rf
      fxp=-1.221774115422d0/rs
      exp=fxp*.75d0
      cpl=g(x)
      cfl=g(y)
      fcp=-cp*cpl
      fcf=-cf*cfl
      ecp=-cp*f(cpl,x)
      ecf=-cf*f(cfl,y)
      cnu=gamma*(ecf-ecp)
      exc=exp+ecp+(fxp+cnu)*ff/gamma
      c1=fxp+cnu
      c2=fcp-cnu
      tcf=(fcf-fcp-4d0*(ecf-ecp)/3d0)*ff
      xu=c1*(1d0+rt)**t+c2+tcf
      xd=c1*(1d0-rt)**t+c2+tcf
c
      else
      xu=0d0
      xd=0d0
      exc=0d0
      endif
c
      u=w(k,ia,1)
      v(k,ia,1)=v(k,ia,1)+xu
      v(k,ia,2)=v(k,ia,2)+xd
      w(k,ia,1)=ro(k,ia,1)*(u+exc)
c     w(k,2)=ro(k,1)*u-3d0*(exc-xu)*(ro(k,1)+ro(k,2))*5d-1
c    &                -3d0*(exc-xd)*(ro(k,1)-ro(k,2))*5d-1
   10 continue
      end
