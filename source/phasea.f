      subroutine phasea(v,tchc,tchs,fcs,mxl,mxlcmp,elvl,tm,ew,ez,ng
     &                 ,dr,xr,meshr,rstr,wk,isr,msg)
c-----------------------------------------------------------------------
c     Construct the Tchebycheff expansion of the phase function.
c     coded by H.Akai, 1983, julich
c     modified by H. Akai, 25 Aug. 1999, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(mx=7)
      real*8 tchc(ng,mxl**2),tchs(ng,mxl**2),v(meshr),dr(meshr)
     &      ,xr(meshr),fcs(3,mxl**2),elvl(ng),tm(ng,ng),rj(mx)
     &      ,dj(mx),rn(mx),dn(mx),rstr(meshr,mxl**2,ng),wk(meshr)
      real*8,allocatable::sn(:,:),cn(:,:),tc1(:),tc2(:),tc(:),td(:)
     &      ,tchr(:,:)
      logical srl,msg,chk
      data c/274.0720442d0/,chk/.false./
      allocate(sn(mx,ng),cn(mx,ng),tc1(ng),tc2(ng),tc(ng),td(ng)
     &        ,tchr(ng,mx))
      srl=isr .eq. 1
      if(mxlcmp .gt. mx) call errtrp(1,'phasea','mxlcmp too large')
      mxj=mxl**2
      lmx=mxlcmp-1
      rtin=xr(meshr-1)
      call clrarr(tchc,ng*mxj)
      call clrarr(tchs,ng*mxj)
      call clrarr(tchr,ng*mx)
      do 20 k=1,ng
      ev=elvl(k)
      if(srl) ev=ev+(ev/c)**2
      call mdbsl(ev*rtin**2,lmx,rj,dj,rn,dn)
      do 20 j=1,mxlcmp
      jj=j
      i=(j-1)**2+1
      call radial(elvl(k),jj,rstr(1,i,k),rxg,dxg,v,v,dr,xr,meshr,isr)
      sn(j,k)=rtin**j*(rtin*rj(j)*dxg-dj(j)*rxg)
      cn(j,k)=rtin**(1-j)*(rtin*rn(j)*dxg-dn(j)*rxg)
      do 20 l=1,ng
   20 tchr(l,j)=tchr(l,j)+rxg*tm(l,k)
      e1=ew-ez
      e2=ew+ez
      call gntcs(e1,ew,ez,tc1,ng)
      call gntcs(e2,ew,ez,tc2,ng)
      do 40 j=1,mxlcmp
      jj=j
      i=(j-1)**2+1
      r1=0d0
      r2=0d0
      do 80 l=1,ng
      r1=r1+tchr(l,j)*tc1(l)
   80 r2=r2+tchr(l,j)*tc2(l)
c     r1=1d0
c     r2=1d0
      e0=e2
c     e0=e1
      if(r1*r2 .lt. 0d0) then
      e0=ew
      de=ez
      do 140 itr=1,4
      call gntcs(e0,ew,ez,tc,ng)
      rr=0d0
      do 150 l=1,ng
  150 rr=rr+tchr(l,j)*tc(l)
      sgn=1d0
      if(rr*r2 .gt. 0d0) sgn=-1d0
      de=5d-1*de
  140 e0=e0+sgn*de
      do 90 itr=1,100
      call gntcds(e0,ew,ez,tc,td,ng)
      rr=0d0
      rd=0d0
      do 100 l=1,ng
      rr=rr+tchr(l,j)*tc(l)
  100 rd=rd+tchr(l,j)*td(l)
      dlt=-rr/rd
      if(abs(dlt) .lt. 1d-5) go to 110
   90 e0=e0+dlt
c
      call errtrp(2,'phasea','e0 not found')
      write(6,'(a,i2)')' l=',j-1
      e0=e2
c
  110 if(e0 .lt. e1 .or. e0 .gt. e2) then
      call errtrp(2,'phasea','illegal e0, converted')
      e0=e2
      endif
      if(msg .and. chk) write(*,'(a,i2,a,f12.7)')' j=',j,'  e0=',e0
      endif
c
      ev=e0
      if(srl) ev=ev+(ev/c)**2
      call mdbsl(ev*rtin**2,jj-1,rj,dj,rn,dn)
      call radial(e0,jj,wk,rxg,dxg,v,v,dr,xr,meshr,isr)
      fcs(1,i)=rtin**(1-j)*(rtin*rn(j)*dxg-dn(j)*rxg)
      fcs(2,i)=rtin**j*(rtin*rj(j)*dxg-dj(j)*rxg)
      fcs(3,i)=e0
      do 40 k=1,ng
      cn(j,k)=(cn(j,k)-fcs(1,i))/(elvl(k)-fcs(3,i))
      sn(j,k)=(sn(j,k)-fcs(2,i))/(elvl(k)-fcs(3,i))
      do 40 l=1,ng
      tchc(l,i)=tchc(l,i)+cn(j,k)*tm(l,k)
   40 tchs(l,i)=tchs(l,i)+sn(j,k)*tm(l,k)
      do 160 j=1,mxlcmp
      i0=(j-1)**2+1
      do 160 i=i0+1,j**2
      do 170 n=1,ng
  170 call equarr(rstr(1,i0,n),rstr(1,i,n),meshr)
      call equarr(tchc(1,i0),tchc(1,i),ng)
      call equarr(tchs(1,i0),tchs(1,i),ng)
  160 call equarr(fcs(1,i0),fcs(1,i),3)
c     write(*,'(1p6e14.6)')
c    & (xr(k),rstr(k,10,1),rstr(k,10,4),rstr(k,10,7),rstr(k,10,11),
c    &  rstr(k,10,15),k=1,meshr-1)
c     stop
      deallocate(sn,cn,tc1,tc2,tc,td,tchr)
      end
