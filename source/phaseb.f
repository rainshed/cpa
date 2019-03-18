      subroutine phaseb(v,tchc,tchs,fcs,mxl,mxlcmp,elvl,tm,ew,ez,ng,dr
     &                 ,xr,meshr,rstr,wk,isr,is,msg)
c-----------------------------------------------------------------------
c     -----------------------------------
c     --- spin-orbit included version ---
c     -----------------------------------
c     Construct the Tchebycheff expansion of the phase function.
c     coded by H.Akai, 1983, julich
c     This version includes the calculation of the spin-orbit coupling
c     constants.
c     coded by H.Akai, 1994, Osaka
c     adopted to spin-orbit version, Nov 1995, Osaka
c     modified by H. Akai, 25 Aug. 1999, Osaka
c     Modified by H. Akai, Tokyo, Jan. 24, 2018.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(mx=7,mx2=mx**2)
      real*8 v(meshr),dr(meshr),xr(meshr),wk(meshr,2*mxl+1)
     &      ,elvl(ng),tm(ng,ng),rj(mx),dj(mx),rn(mx),dn(mx)
     &      ,fcs(3,mxl**2),tchc(ng,mxl**2)
     &      ,tchs(ng,mxl**2),rstr(meshr,mxl**2,ng)
      real*8,allocatable::sn(:,:),cn(:,:),tchr(:,:)
     &      ,tc1(:),tc2(:),tc(:),td(:)
      logical srl,msg,chk,test
      data c/274.0720442d0/,chk/.false./,test/.false./
      allocate(sn(mx2,ng),cn(mx2,ng),tchr(ng,mx2),tc1(ng),tc2(ng),tc(ng)
     &        ,td(ng))
c     if(i .eq. 2) then
c     test=.true.
c     else
c     test=.false.
c     endif
      if(mxlcmp .gt. mx) call errtrp(1,'phaseb','mxlcmp too large')
      if(test) call errtrp(3,'phaseb','spin-orbit suppressed')
      srl=isr .eq. 1
      spn=(-1d0)**(is-1)
      mxj=mxl**2
      lmax=mxlcmp-1
      mxjcmp=mxlcmp**2
      rtin=xr(meshr-1)
      v0=v(meshr)-v(meshr-1)
      if(chk) write(*,*)
      if(test) v0=0d0
      call clrarr(tchc,ng*mxj)
      call clrarr(tchs,ng*mxj)
      call clrarr(tchr,ng*mx2)
      do 102 k=1,meshr-1
  102 wk(k,1)=v(k)*xr(k)
      call diffn(wk(1,1),wk(1,2),xr,dr,meshr-1)
      do 100 k=1,meshr-1
  100 wk(k,2)=(wk(k,2)-v(k))/xr(k)**2
      do 20 k=1,ng
      ev=elvl(k)
      if(srl) ev=ev+(ev/c)**2
      call mdbsl(ev*rtin**2,lmax,rj,dj,rn,dn)
      do 10 m=-lmax,lmax
      mm=3+lmax+m
      wk(meshr,mm)=v(meshr)
      fct=spn*dble(m)/c**2
c     --- spin moment 1/2 and factro 2 appearing in spin-orbit
c         coupling cancels to each other, remaining no extra
c         factor.
c     fct=spn*dble(-m)/c**2
c     write(*,*)k,m,fct
      if(test) fct=0d0
c     if(abs(m) .eq. 2) fct=0d0
      if(srl) then
      do 110 kk=1,meshr-1
  110 wk(kk,mm)=v(kk)+fct*wk(kk,2)/(1d0-(v(kk)-elvl(k))/c**2)**2
      else
      do 130 kk=1,meshr-1
  130 wk(kk,mm)=v(kk)+fct*wk(kk,2)
      endif
   10 continue
      do 20 j=1,mxjcmp
      jj=sqrt(dble(j)-5d-1)+1
      m=j-1-jj*(jj-1)
      mm=3+lmax+m
      call radial(elvl(k),jj,rstr(1,j,k),rxg,dxg,v,wk(1,mm),dr,xr
     &           ,meshr,isr)
      if(.false.)then
      if(srl) then
      do 120 kk=1,meshr-1
  120 wk(kk,1)=wk(kk,2)*(rstr(kk,j,k)/xr(kk))**2
     &        /(1d0-(v(kk)-elvl(k))/c**2)**2
      else
      do 122 kk=1,meshr-1
  122 wk(kk,1)=wk(kk,2)*(rstr(kk,j,k)/xr(kk))**2
      endif
c     sum=0d0
c     do 124 kk=1,meshr-3,2
c     sum=sum+2d0*(wk(kk,1)*xr(kk)**2*dr(kk)
c    &       +4d0*wk(kk+1,1)*xr(kk+1)**2*dr(kk+1)
c    &       +wk(kk+2,1)*xr(kk+2)**2*dr(kk+2))/3d0/c**2
c 124 write(*,'(3x,i4,1p5e15.6)')kk+1,xr(kk+1),rstr(kk+1,j,k)
c    &    ,wk(kk+1,2)*xr(KK+1)**2,wk(kk+1,mm)*xr(kk+1),sum
      slprm=2d0*fintgr(wk,dr,xr,meshr)/c**2
c     slprm is increasing function with respect to energy because
c     the nodes of higher energy radial wave function come
c     more inside.
      write(*,'(1x,a,f10.5,2i3,f10.5)')'  e,l,m,sl=',
     &     elvl(k),jj,m,slprm
      endif
      dxg=dxg+spn*dble(m)*v0*rxg/c**2
      sn(j,k)=rtin**jj*(rtin*rj(jj)*dxg-dj(jj)*rxg)
      cn(j,k)=rtin**(1-jj)*(rtin*rn(jj)*dxg-dn(jj)*rxg)
      do 20 l=1,ng
   20 tchr(l,j)=tchr(l,j)+rxg*tm(l,k)
      e1=ew-ez
      e2=ew+ez
      call gntcs(e1,ew,ez,tc1,ng)
      call gntcs(e2,ew,ez,tc2,ng)
      do 30 j=1,mxjcmp
      jj=sqrt(dble(j)-5d-1)+1
      m=j-1-jj*(jj-1)
      mm=3+lmax+m
      r1=0d0
      r2=0d0
      do 40 l=1,ng
      r1=r1+tchr(l,j)*tc1(l)
   40 r2=r2+tchr(l,j)*tc2(l)
      e0=e2
      if(r1*r2 .lt. 0d0) then
      e0=ew
      de=ez
      do 50 itr=1,4
      call gntcs(e0,ew,ez,tc,ng)
      rr=0d0
      do 60 l=1,ng
   60 rr=rr+tchr(l,j)*tc(l)
      sgn=1d0
      if(rr*r2 .gt. 0d0) sgn=-1d0
      de=5d-1*de
   50 e0=e0+sgn*de
      do 70 itr=1,100
      call gntcds(e0,ew,ez,tc,td,ng)
      rr=0d0
      rd=0d0
      do 80 l=1,ng
      rr=rr+tchr(l,j)*tc(l)
   80 rd=rd+tchr(l,j)*td(l)
      dlt=-rr/rd
      if(abs(dlt) .lt. 1d-5) go to 90
   70 e0=e0+dlt
c
      call errtrp(2,'phaseb','e0 not found')
      write(*,'(a,i2)') '      for j=',j
      e0=e2
c
   90 if(e0 .lt. e1 .or. e0 .gt. e2) then
      call errtrp(2,'phaseb','illegal e0, converted')
      e0=e2
      endif
      if(msg .and. chk) then
      call errtrp(3,'phaseb',' ')
      write(*,'(6x,a,i1,a,f12.7)')'j=',j,'  e0=',e0
      endif
      endif
c
      ev=e0
      if(srl) then
      ev=ev+(ev/c)**2
      wk(meshr,mm)=v(meshr)
      fct=spn*dble(m)/c**2
      if(test) fct=0d0
c     if(abs(m) .eq. 2) fct=0d0
      do 140 kk=1,meshr-1
  140 wk(kk,mm)=v(kk)+fct*wk(kk,2)/(1d0-(v(kk)-e0)/c**2)**2
      endif
      call mdbsl(ev*rtin**2,jj-1,rj,dj,rn,dn)
      call radial(e0,jj,wk,rxg,dxg,v,wk(1,mm),dr,xr,meshr,isr)
      dxg=dxg+spn*dble(m)*v0*rxg/c**2
      fcs(1,j)=rtin**(1-jj)*(rtin*rn(jj)*dxg-dn(jj)*rxg)
      fcs(2,j)=rtin**jj*(rtin*rj(jj)*dxg-dj(jj)*rxg)
      fcs(3,j)=e0
      if(chk)then
      if(srl) then
      do 150 kk=1,meshr-1
  150 wk(kk,mm)=wk(kk,2)*(wk(kk,1)/xr(kk))**2
     &        /(1d0-(v(kk)-e0)/c**2)**2
      else
      do 160 kk=1,meshr-1
  160 wk(kk,mm)=wk(kk,2)*(wk(kk,1)/xr(kk))**2
      endif
      slprm=2d0*fintgr(wk(1,mm),dr,xr,meshr)/c**2
      write(*,'(1x,a,f10.5,2i3,f10.5)')'  e,l,m,sl=',
     &     e0,jj-1,m,slprm
      endif
      do 30 k=1,ng
      cn(j,k)=(cn(j,k)-fcs(1,j))/(elvl(k)-fcs(3,j))
      sn(j,k)=(sn(j,k)-fcs(2,j))/(elvl(k)-fcs(3,j))
      do 30 l=1,ng
      tchc(l,j)=tchc(l,j)+cn(j,k)*tm(l,k)
   30 tchs(l,j)=tchs(l,j)+sn(j,k)*tm(l,k)
c     write(*,*)'phaseb'
c     write(*,'(1x,a,16f10.5)')' fcs(1)=',(fcs(1,j),j=1,mxlcmp**2)
c     write(*,'(1x,a,16f10.5)')' fcs(2)=',(fcs(2,j),j=1,mxlcmp**2)
c     write(*,'(1x,a,16f10.5)')' fcs(3)=',(fcs(3,j),j=1,mxlcmp**2)
c     write(*,'(1x,a,16f10.5)')' fcs(4)=',(fcs(4,j),j=1,mxlcmp**2)
c     write(*,'(1x,a,(5f10.5))')' cn=',(cn(5,k),k=1,ng)
c     write(*,'(1x,a,(5f10.5))')' sn=',(sn(5,k),k=1,ng)
c     write(*,'(1x,a,(5f10.5))')' tchc=',(tchc(k,5),k=1,ng)
c     write(*,'(1x,a,(5f10.5))')' tchc=',(tchs(k,5),k=1,ng)
      deallocate(sn,cn,tchr,tc1,tc2,tc,td)
      end
