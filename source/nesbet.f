      subroutine nesbet(t,tcpa,f,cj,conc,dmpc,ne,ncmp,ntyp,ncmpx,natm
     &                 ,itype,convrg,mxl,alloy,lmxtyp)
c-----------------------------------------------------------------------
c     Single step iteration of coherent t-matrix by Nesbet's algorithm.
c     Ref. R.K. Nesbet, Phys. Rev. B45 (1992) 13234.
c     f; sum (1/(-1/t + b(struct) +i sqrt(e)))
c     The self consistent equation of CPA is
c     c_A/((-1/t_A)-J)+c_B/((-1/t_B)-J)=1/((-1/t)-J).
c     Nesbet's algorithm is the following.
c     define m'=(-1/tcpa)_new, m=(-1/tcpa)_old, m_i=(-1/t_i),
c     and T=Sum_k [(-1/tcpa)_old - g_struct]^(-1).
c     The iteration process is
c     1/m'=1/m + Sum_i (1/m) c_i*[1/(m-m_i)-T]^(-1) (1/m)
c     Coded by H.Akai, 14 Sep 1996, Osaka
c     Spin-orbit version, 18 April 1997, Osaka
c     Full matrix version, 13 August 1997, Duisburg
c     Revised so as to avoid divergence occuring at m-m_i=0.
c                 by H. Akai, 23 Jan. 2005, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(mxindx=36)
      real*8 conc(ncmpx)
      complex*16 t(mxl**2,ncmpx,ne),tcpa(mxl**2,mxl**2,natm,ne)
     &          ,f(mxl**2,mxl**2,natm,ne),cj(mxl**2,mxl**2,natm,ne)
      complex*16,allocatable::wk(:,:,:)
      integer ncmp(ntyp),itype(natm),lmxtyp(ntyp),indx(mxindx)
      logical convrg(ne),alloy
      data tol/1d-8/,zero/1d-8/
      allocate(wk(mxl**2,mxl**2,3))
      if(mxl**2 .gt. mxindx) call errtrp(1,'nesbet','mxindx too small')
      mmxl=mxl**2
      mf=mmxl**2
c
c     --- calculate a new coherent t-matrix --
c     --- modify the old one by use of the new one
c
c     call dspcpa(tcpa(1,1,ne),mmxl)
c     write(*,'(1x,a)')'diagonal...nesbet'
      do 10 k=1,ne
c     write(*,'(1x,a,i2,a,l1)')'k=',k,' convrg=',convrg(k)
      if(.not. convrg(k)) then
      error=0d0
c     --- cj=t-1/f is calculated (f=T).
      do 20 i=1,natm
      mmx=(lmxtyp(itype(i))+1)**2
      do 30 mc=1,mmx
      do 30 mr=1,mmx
   30 wk(mr,mc,1)=f(mr,mc,i,k)
      call cinvrn(wk,cj(1,1,i,k),mmxl,mmx,1)
c     --- cj=m-1/T, where m=1/t
      do 20 mc=1,mmx
      do 20 mr=1,mmx
   20 cj(mr,mc,i,k)=tcpa(mr,mc,i,k)-cj(mr,mc,i,k)
      if(alloy) then
      do 40 i=1,natm
      n=itype(i)
      mmx=(lmxtyp(n)+1)**2
      do 90 j=1,ncmp(n)
      call jip(n,j,jn)
c     --- if any of concentrations is 1, that site is pure.
      if(abs(conc(jn)-1d0) .lt. zero) go to 40
   90 continue
      call clrarc(wk(1,1,1),mf)
      do 50 j=1,ncmp(n)
      call jip(n,j,jn)
c     --- wk2=m-m_i
      do 60 mc=1,mmx
      do 70 mr=1,mmx
   70 wk(mr,mc,2)=tcpa(mr,mc,i,k)
   60 wk(mc,mc,2)=wk(mc,mc,2)-t(mc,jn,k)
c     --- wk3=-(m-m_i)*T
      do 80 mc=1,mmx
      do 80 mr=1,mmx
      wk(mr,mc,3)=(0d0,0d0)
      do 80 l=1,mmx
   80 wk(mr,mc,3)=wk(mr,mc,3)-wk(mr,l,2)*f(l,mc,i,k)
c     --- wk3=1-(m-m_i)T    
      do 150 mc=1,mmx
  150 wk(mc,mc,3)=wk(mc,mc,3)+1d0   
c     --- wk2=(1/(1-(m-m_i)T))*(m-m_i)
      call ludcmp(wk(1,1,3),mmx,mmx,indx,parity)
      do 160 mc=1,mmx
  160 call lubksb(wk(1,1,3),mmx,mmx,indx,wk(1,mc,2))
c     --- wk1=sum(c_i*(m-m_i)/(1-T*(m-m_i)))
      do 50 mc=1,mmx
      do 50 mr=1,mmx
   50 wk(mr,mc,1)=wk(mr,mc,1)+conc(jn)*wk(mr,mc,2)
c     --- wk3=1/m
      do 120 mc=1,mmx
      do 120 mr=1,mmx
  120 wk(mr,mc,2)=tcpa(mr,mc,i,k)
      call cinvrn(wk(1,1,2),wk(1,1,3),mmxl,mmx,1)
c     --- wk3=1/m
c     --- wk2=sum(...) * (1/m)
      do 130 mc=1,mmx
      do 130 mr=1,mmx
      wk(mr,mc,2)=(0d0,0d0)
      do 130 l=1,mmx
  130 wk(mr,mc,2)=wk(mr,mc,2)+wk(mr,l,1)*wk(l,mc,3)
c     --- wk1=(1/m) * sum(...) * (1/m)
      do 140 mc=1,mmx
      do 140 mr=1,mmx
      wk(mr,mc,1)=(0d0,0d0)
      do 140 l=1,mmx
  140 wk(mr,mc,1)=wk(mr,mc,1)+wk(mr,l,3)*wk(l,mc,2)
c     --- check the maximum correction to be made
      chk=0d0
      do 100 mc=1,mmx
      do 100 mr=1,mmx
      sn=abs(wk(mr,mc,3))
  100 chk=max(chk,abs(wk(mr,mc,1))*sn/(1d0+sn**2))
      if(chk .gt. tol) then
c     --- too big a correction make iteration diverging
c     --- reduce dmpc to avoid divergence 
c     --- following procedure is not final
      dmp=dmpc/max(chk,1d0)
      do 110 mc=1,mmx
      do 110 mr=1,mmx
  110 wk(mr,mc,3)=wk(mr,mc,3)+dmp*wk(mr,mc,1)
      call cinvrn(wk(1,1,3),tcpa(1,1,i,k),mmxl,mmx,1)
      endif
      error=max(chk,error)
   40 continue
      endif
      if(error .lt. tol) convrg(k)=.true. 
      endif
   10 continue
c     call dspcpa(tcpa(1,1,ne),mmxl)
c     end
c     subroutine dspcpa(t,m)
c     complex*16 t(m,m)
c     write(*,'((1x,9f12.7))')(abs(t(j,j)),j=1,m)
      deallocate(wk)
      end
