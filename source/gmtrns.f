      subroutine gmtrns(gm,nk,mxla,mxlb,msiz,iap,ibp,proc,uu,mxl)
c----------------------------------------------------------------------
c     Given a matrixi gm, this program transform it either from
c     spherical harmonics to real harmonics representation (proc='s2r')
c     or the other way (proc='r2s'). msiz is the dimension of gm.
c     Only the nda by ndb sub-block starting from gm(iap,ibp)
c     is transformed. 
c     Originally coded by Manabu Takahashi, IMR.,Nov,1994
c     Modified by H. Akai, 22 Aug. 1997, Duisburg
c     Completely modified by H. Akai, 25 Aug. 1999, Osaka
c     Newly Coded by H. Akai, 9 April 2014, Tokyo
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 gm(msiz,msiz,nk),uu(2*mxl-1,2*mxl-1,2)
      complex*16,allocatable::wk(:,:)
      character proc*(*)
      data zero/1d-10/
      nl=max(mxla,mxlb)**2
      allocate(wk(nl,nl),stat=ierr)
      if(ierr .ne. 0) call errtrp(1,'gmtrns','allocaltion fails')
      ia=iap-1
      ib=ibp-1
      nda=mxla**2
      ndb=mxlb**2
      if(proc .eq. 'S2R' .or. proc .eq. 's2r') then
      il=1
      else if(proc .eq. 'R2S' .or. proc .eq. 'r2s') then
      il=2
      else
      call errtrp(1,'gmtrns','illegal proc')
      endif
      ir=3-il
      do 10 kk=1,nk
      call clrarc(wk,nl**2)
      do 20 lb=0,mxlb-1
      mb=lb*(lb+1)+1
      mu=mxl-mb
      do 20 j=mb-lb,mb+lb
      do 20 i=mb-lb,mb+lb
      if(abs(uu(i+mu,j+mu,ir)) .gt. zero) then
      do 30 k=1,nda
   30 wk(k,j)=wk(k,j)+gm(k+ia,i+ib,kk)*uu(i+mu,j+mu,ir)
      endif
   20 continue
      do 40 j=1,ndb
      do 40 i=1,nda
   40 gm(i+ia,j+ib,kk)=(0d0,0d0)
      do 10 la=0,mxla-1
      ma=la*(la+1)+1
      mu=mxl-ma
      do 10 j=ma-la,ma+la
      do 10 i=ma-la,ma+la
      if(abs(uu(i+mu,j+mu,il)) .gt. zero) then
      do 50 k=1,ndb
   50 gm(i+ia,k+ib,kk)=gm(i+ia,k+ib,kk)+uu(i+mu,j+mu,il)*wk(j,k)
      endif
   10 continue
      deallocate(wk,stat=ierr)
      end
