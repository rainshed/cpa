      subroutine green(g,pf,phase,str,ck,f,cj,t,mxl,ntyp,ncmp,iatm
     &                 ,kmx,natm,ncmpx,lmxtyp)
c----------------------------------------------------------------------
c     This program calculate the componet Green fucntion together
c     with their phase.
c     Coded by H. Akai 14 August, 1997, Duisburg
c     Modified by H. Akai, 25 Aug. 1999, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 phase(mxl**2,ncmpx,kmx)
     &          ,str(mxl**2,ncmpx,kmx),cj(mxl**4,natm,kmx)
     &          ,t(mxl**2,ncmpx,kmx),f(mxl**4,natm,kmx)
     &          ,ck(mxl**4,natm,kmx),g(mxl**2,ncmpx,kmx)
     &          ,pf(mxl**2,mxl**2,ncmpx,kmx)
      complex*16,allocatable::wk(:,:)
      integer iatm(ntyp),ncmp(ntyp),lmxtyp(ntyp)
      call sbtime(4,0)
      allocate(wk(mxl**4,3))
      mmxl=mxl**2
      mf=mmxl**2
      do 10 k=1,kmx
c     write(*,'((1x,i3,9f12.6))')k,(imag(str(i,1,k)),i=1,9)
c     write(*,'((1x,i3,9f12.6))')k,(imag(ck(mmxl*(i-1)+i,1,k)),i=1,9)
      do 10 i=1,ntyp
      mmx=(lmxtyp(i)+1)**2
      n=iatm(i)
c     do 20 m=1,mf
c  20 wk(m,2)=f(m,n,k)
c     do 30 m=1,mmx
c     mm=mmxl*(m-1)+m
c  30 wk(mm,2)=log(wk(mm,2))+log(-ck(mm,n,k))
      do 10 j=1,ncmp(i)
      call jip(i,j,ji)
      do 40 m=1,mf
   40 wk(m,1)=-cj(m,n,k)
      do 50 m=1,mmx
      mm=mmxl*(m-1)+m
   50 wk(mm,1)=wk(mm,1)+t(m,ji,k)
      call prdmtc(wk(1,1),f(1,n,k),wk(1,3),mmxl,mmx)
      call cinvrn(wk(1,3),wk(1,1),mmxl,mmx,1)
      call prdmtc(f(1,n,k),wk(1,1),wk(1,2),mmxl,mmx)
      do 10 m=1,mmx
      mm=mmxl*(m-1)+m
      phase(m,ji,k)=str(m,ji,k)+log(-ck(mm,n,k))+log(wk(mm,3))
c     phase(m,ji,k)=str(m,ji,k)
c     phase(m,ji,k)=log(-ck(mm,n,k))+log(wk(mm,3))
c     phase(m,ji,k)=str(m,ji,k)+log(wk(mm,3))
   10 g(m,ji,k)=g(m,ji,k)+t(m,ji,k)*(1d0-t(m,ji,k)
     &           *wk(mm,2))*pf(m,m,ji,k)
      deallocate(wk)
      call sbtime(4,1)
      end
