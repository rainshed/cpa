      subroutine cryrho(ro,dr,xr,msr,a,rox,xrx,meshx,nei,dist
     &                 ,ityp,deg,wk,icmp,ncmp,ntyp,ncmpx)
c-----------------------------------------------------------------------
c     -----------------------
c     --- KKR-CPA version ---
c     -----------------------
c     Construct charge density after Mattheiss' prescription.
c     coded by H.Akai, Jan. 1986, Osaka
c     KKR-CPA implemented by H.Akai, 7 Sep. 1996, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 ro(msr),dr(msr),xr(msr),rox(meshx,ncmpx)
     &      ,xrx(meshx,ncmpx),dist(nei),deg(nei),wk(msr)
      integer ityp(nei),ncmp(ntyp)
      if(ityp(1) .eq. 0) call errtrp(1,'cryrho','no type assigned')
c
      kmx=(msr-1)/2
      wk(1)=xr(1)
      wk(kmx+1)=0d0
      do 10 k=2,kmx
      n=2*k-3
      wk(k)=xr(n+2)
   10 wk(kmx+k)=wk(kmx+k-1)+(dr(n)*xr(n)*ro(n)
     &      +4d0*dr(n+1)*xr(n+1)*ro(n+1)+dr(n+2)*xr(n+2)*ro(n+2))/3d0
      i=ityp(1)
      call jip(i,icmp,icmpi)
      do 20 k=1,meshx
   20 rox(k,icmpi)=rox(k,icmpi)+polint(xrx(k,icmpi),xr,ro,msr)
      do 30 n=2,nei
      i=ityp(n)
      if(i .eq. 0) return
      do 30 ic=1,ncmp(i)
      call jip(i,ic,ici)
      d=a*dist(n)
      factor=deg(n)*5d-1/d
      do 30 k=1,meshx
      x1=d+xrx(k,ici)
      x2=abs(d-xrx(k,ici))
      sum=0d0
      if(x1 .lt. xr(msr-2)) sum=sum+polint(x1,wk(1),wk(kmx+1),kmx)
      if(x2 .lt. xr(msr-2)) sum=sum-polint(x2,wk(1),wk(kmx+1),kmx)
      rox(k,ici)=rox(k,ici)+factor*sum/xrx(k,ici)
   30 continue
      return
      end
