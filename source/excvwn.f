      subroutine excvwn(ro,v,w,dr,xr,msr,mna,meshr)
c----------------------------------------------------------------------
c     +----------------------------------+
c          SDF parametrization by VWN
c     +----------------------------------+
c     'uxcor' coded by  M.Mannien is called in this program
c     coded by H.Akai, 1984, Juelich
c     very minor revision by H.Akai, 23 Nov. 1995, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 ro(msr,mna,2),v(msr,mna,2),w(msr,mna,*)
     &      ,dr(msr,mna),xr(msr,mna)
      do 10 ia=1,mna
      do 10 k=1,meshr
c
      if(ro(k,ia,1) .gt. 1d-20 .and.
     &       .not. abs(ro(k,ia,2)) .gt. ro(k,ia,1)) then
      rs=(3d0/ro(k,ia,1))**(1d0/3d0)
      rt=ro(k,ia,2)/ro(k,ia,1)
      rt=min(1d0,max(-1d0,rt))
      call uxcor(rs,rt,xu,xd,exc)
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
