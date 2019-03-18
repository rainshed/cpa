      subroutine cpaitr(t,tcpa,f,cj,conc,dmpc,ne,ncmp,ntyp,ncmpx,natm
     &                 ,itype,convrg,mxl,alloy,lmxtyp)
c-----------------------------------------------------------------------
c     Single step iteration of coherent t-matrix
c     f; sum (1/(-1/t + b(struct) +i sqrt(e)))
c     The self consistent equation of CPA is
c     c_A/((-1/t_A)-J)+c_B/((-1/t_B)-J)=1/((-1/t)-J).
c     Coded by H.Akai, 14 Sep 1996, Osaka
c     Spin-orbit version, 18 April 1997, Osaka
c     Full matrix version, 13 August 1997, Duisburg
c     Modified by H. Akai, 25 Aug. 1999, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 conc(ncmpx)
      complex*16 t(mxl**2,ncmpx,ne),tcpa(mxl**2,mxl**2,natm,ne)
     &          ,f(mxl**2,mxl**2,natm,ne),cj(mxl**2,mxl**2,natm,ne)
      complex*16,allocatable::wk(:,:,:)
      integer ncmp(ntyp),itype(natm),lmxtyp(ntyp)
      logical convrg(ne),alloy,notyet
      character mark*201
      data tol/1d-6/, zero/1d-8/
      allocate(wk(mxl**2,mxl**2,3))
      mmxl=mxl**2
      mf=mmxl**2
c
c     --- calculate a new coherent t-matrix --
c     --- modify the old one by use of the new one
c
c     call dspcpa(tcpa(1,1,ne),mmxl)
c     write(*,*)'diagonal...cpaitr'
      do 10 k=1,ne
      if(.not. convrg(k)) then
      error=0d0
c     --- cj=t-1/f is calculated.
c     if(k .eq. ne)
c    &  write(*,'(/(1x,9f10.5))')
c    &  ((abs(f(mr,mc,1,k)),mr=1,mmxl),mc=1,mmxl)
      do 20 i=1,natm
      mmx=(lmxtyp(itype(i))+1)**2
      do 22 mc=1,mmx
      do 22 mr=1,mmx
   22 wk(mr,mc,1)=f(mr,mc,i,k)
      call cinvrn(wk,cj(1,1,i,k),mmxl,mmx,1)
      do 20 mc=1,mmx
      do 20 mr=1,mmx
   20 cj(mr,mc,i,k)=tcpa(mr,mc,i,k)-cj(mr,mc,i,k)
c     if(k .eq. ne)
c    &  write(*,'(/(1x,9f10.5))')
c    &     (abs(cj(mr,mc,1,k)),mr=1,mmxl),mc=1,mmxl)
      if(alloy) then
      do 30 i=1,natm
      n=itype(i)
      mmx=(lmxtyp(n)+1)**2
      do 90 j=1,ncmp(n)
      call jip(n,j,jn)
c     --- if any of concentrations is 1, that site is pure.
      if(abs(conc(jn)-1d0) .lt. zero) go to 30
   90 continue
      call clrarc(wk,mf)
      do 40 j=1,ncmp(n)
      call jip(n,j,jn)
      do 50 mc=1,mmx
      do 50 mr=1,mmx
   50 wk(mr,mc,2)=-cj(mr,mc,i,k)
      do 60 mr=1,mmx
   60 wk(mr,mr,2)=wk(mr,mr,2)+t(mr,jn,k)
      call cinvrn(wk(1,1,2),wk(1,1,3),mmxl,mmx,1)
      do 40 mc=1,mmx
      do 40 mr=1,mmx
   40 wk(mr,mc,1)=wk(mr,mc,1)+conc(jn)*wk(mr,mc,3)
      call cinvrn(wk(1,1,1),wk(1,1,3),mmxl,mmx,1)
      chk=0d0
      do 70 mc=1,mmx
      do 70 mr=1,mmx
      wk(mr,mc,1)=cj(mr,mc,i,k)+wk(mr,mc,3)-tcpa(mr,mc,i,k)
      sn=abs(tcpa(mr,mc,i,k))
      err=abs(wk(mr,mc,1))*sn/(1d0+sn**2)
c     if(sn .lt. 1d-6) then
c     err=sn
c     else
c     err=abs(wk(m,1))/sn
c     endif
      if(err .gt. chk) then
      mcer=mc
      mrer=mr
      chk=err
      endif
   70 continue
c  70 chk=max(chk,abs(wk(m,1))*sn/(1d0+sn**2))
c  70 chk=max(chk,abs(wk(m,1))/max(1d0,abs(tcpa(m,i,k))))
c     if(k .eq. ne .and. i .eq. 1)
c    &  write(*,'(/1x,9f10.5)')(abs(tcpa(m,1,k)),m=1,mf,10)
c    &   write(*,'(/1x,i3,3f12.7)')mer,chk,tcpa(mer,i,k)
      if(chk .gt. tol) then
      do 80 mc=1,mmx
      do 80 mr=1,mmx
   80 tcpa(mr,mc,i,k)=tcpa(mr,mc,i,k)+dmpc*wk(mr,mc,1)
c     if(k .eq. ne .and. i .eq. 1)
c    &  write(*,'(1x,9f10.5)')(abs(tcpa(m,1,k)),m=1,mf,10)
c    &  write(*,'(1x,9f10.5)')(abs(wk(m,1)),m=1,mf,10)
c    &  write(*,'((1x,9f10.5))')(abs(tcpa(m,1,k)),m=1,mf)
c     if(k .eq. ne .and. i .eq. 1)
c     if(k .eq. ne)
c    &   write(*,'(1x,3i3,3f12.7)')i,mrer,mcer,chk,tcpa(mrer,mcer,i,k)
      endif
      error=max(chk,error)
   30 continue
      endif
      if(error .lt. tol) convrg(k)=.true.
      endif
   10 continue
c     write(*,'(1x,35l1)')convrg
      mark=' '
      notyet=.false.
      do 100 k=1,ne
      if(.not. convrg(k)) then
      mark(k:k)='*'
      notyet=.true.
      endif
  100 continue
c     if(notyet) write(*,'(1x,a)')mark(1:ne)
      deallocate(wk)
      end
