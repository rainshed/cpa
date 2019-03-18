      subroutine spckkr(g,ff,e,kmx,cm,tchc,tchs,fcs,ew,ez
     &                 ,mxl,ng,iblk,natm,ntyp,itype,vc
     &                 ,isr,a,nk,np,tch,pexf,prr,hh,ncub,clks
     &                 ,last,detl,wtkp,nd,t,pf,str,tc,iwtyp
     &                 ,ids,length,gfree,ess,ncmp,ncmpx,conc
     &                 ,tcpa,phase,korder,convrg,cnvq,urotat,uu
     &                 ,irotat,isymop,ck,cj,iatm,ls,spctrl,nk3
     &                 ,lmxtyp,mxlcmp,msiz,lmxblk,iatmp,itblk,its)
c-----------------------------------------------------------------------
c     ------------------------------------
c     --- KKR-CPA + spin-orbit version ---
c     ------------------------------------
c     Calculates all needed in KKR
c     coded by H.Akai, Aug. 1986, Juelich
c     latest revision by H.Akai, Feb 1996, Osaka
c     KKR CPA implemented by H.Akai, 14 Sep 1996, Osaka
c     spin-orbit version, 18 Apri 1997, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
c     --- Note: length is set in the main program "specx".
      complex*16 g(mxl**2,ncmpx,kmx),e(kmx),detl(kmx)
     &          ,tch((2*mxl-1)**2,ng,nk,nd)
     &          ,pf(mxl**2,mxl**2,ncmpx,kmx)
     &          ,t(mxl**2,ncmpx,kmx),pexf(np,nk,nd),tc(ng,kmx)
     &          ,str(mxl**2,ncmpx,kmx),ff(mxl**2,mxl**2,natm,kmx)
     &          ,gfree(kmx),ess(kmx,(2*mxl-1)**2),cunit
c      --- off diagonal elements are also considered in this version.
     &          ,tcpa(mxl**2,mxl**2,natm,kmx)
     &          ,phase(mxl**2,ncmpx,kmx)
     &          ,ck(mxl**2,mxl**2,natm,kmx),cj(mxl**2,mxl**2,natm,kmx)
     &          ,urotat((2*mxl-1)**2,mxl,24),spctrl(kmx,nk3)
     &          ,uu(2*mxl-1,2*mxl-1,2)
      complex*16,allocatable::gi(:),gm(:),det(:)
c
      real*8     tchc(ng,mxl**2,ncmpx),tchs(ng,mxl**2,ncmpx)
     &          ,cm(ng,mxl**2,ncmpx)
     &          ,prr(np,nk),hh(np,(2*mxl-1)**2,nk)
     &          ,fcs(3,mxl**2,ncmpx),clks(last)
     &          ,wtkp(nk),conc(ncmpx)
c
      integer    itype(natm),iblk(natm,natm),ncub(last)
     &          ,iprint(10),iwtyp(ntyp),lmxblk(nd),itblk(5,its)
     &          ,ncmp(ntyp),korder(nk),irotat(natm,48),isymop(48)
     &          ,iatm(ntyp),lmxtyp(ntyp),mxlcmp(ncmpx),iatmp(natm)
c
      logical    convrg(kmx),start,ptpure,pure,alloy,cnv,poor,nesitr
      logical,allocatable::convr2(:)
      character*1  lsymbl(6)
      data iprint/0,0,0,0,0,0,0,0,0,0/
     &    ,c/274.0720442d0/, cunit/(0d0,1d0)/
     &    ,dmpc0/1d0/, itrmx/50/, nesitr/.false./
c    &    ,dmpc0/2d-1/, itrmx/50/, nesitr/.false./
     &    ,lsymbl/'s','p','d','f','g','h'/
c
c        --- print switch ---
c      iprint
c        1     negative inverse t-matrix
c        2     not used
c        3     log(det green's function)
c        4     print s and c function
c        5     not used
c        6     not used
c        7     not used
c        8     not used
c        9     not used
c       10     not used
c
      allocate(convr2(kmx))
      pi=4d0*atan(1d0)
      nk1=nk-nk3
      unit=2d0*pi/a
      mmxl=mxl**2
      if(isr .eq. 1) then
      relfct=1d0/c**2
      else
      relfct=0d0
      endif
      start=cnvq .gt. 10d0
c     poor=cnvq .gt. 1d0
      poor=cnvq .gt. 1d-2
c     poor=cnvq .gt. 1d-4
c     poor=cnvq .gt. 0d0
c
      if(iprint(4) .ge. 1) then
      do 10 i=1,ncmpx
c     do 10 i=4,4
   10 call prntcs(tchc(1,1,i),tchs(1,1,i),fcs(1,1,i)
     &           ,ew,ez,ng,mxlcmp(i),201,isr)
      endif
c
      do 20 i=1,ncmpx
      do 30 k=1,kmx
      call cgntcs(e(k),ew,ez,tc(1,k),ng)
   30 call cpshft(e(k),tchc(1,1,i),tchs(1,1,i),cm(1,1,i),fcs(1,1,i)
     &          ,g(1,i,k),tc(1,k),mxl,mxlcmp(i),ng,t(1,i,k)
     &          ,pf(1,1,i,k),str(1,i,k),isr,ids)
c
c     --- offset gs(k) by real linear function in order to
c         keep the precision of the complex energy integration.
      do 32 l=1,mxlcmp(i)**2
      a1=(dble(g(l,i,kmx))-dble(g(l,i,1)))
     &   /(dble(e(kmx))-dble(e(1)))
      a0=dble(g(l,i,1))
      do 32 k=1,kmx
   32 g(l,i,k)=g(l,i,k)-a1*(e(k)-dble(e(1)))-a0
c     call clrarc(g,mmxl*ncmpx*kmx)
c     --- print out negative inverse t-matrix ---
      if(iprint(1) .ge. 1) then
      mxdsp=min(5,mxlcmp(i))
      write(*,1100)i,(lsymbl(l),l=1,mxdsp)
 1100 format(/'   -1/t for component',i2/t4,'k',t18,a1,4(20x,a1)) 
      write(*,'(1x,6(''-''),5a)')('---------------------',l=1,mxdsp)
      do 33 k=1,kmx
   33 write(*,'(1x,i3,3x,5(1p,2e10.2,1x))')
     &                  k,(t(l*(l-1)+1,i,k),l=1,mxdsp)
      end if
   20 continue
c
c     --- Starting values for tcpa is the average t-matrix.
c     --- tcpa is initialized only when "start" is true.
      pure=.true.
      do 80 n=1,natm
      i=itype(n)
      ptpure=.false.
      do 82 j=1,ncmp(i)
      call jip(i,j,ji)
   82 if(abs(conc(ji)-1d0) .lt. 1d-10) ptpure=.true.
      pure=pure .and. ptpure
      if(start .or. poor .or. ptpure) then
      mmx=(lmxtyp(i)+1)**2
      do 86 l1=1,mmx
      do 86 l2=1,mmx
      do 86 k=1,kmx
   86 tcpa(l1,l2,n,k)=(0d0,0d0)
      do 90 j=1,ncmp(i)
      call jip(i,j,ji)
      do 90 l=1,mmx
      do 90 k=1,kmx
   90 tcpa(l,l,n,k)=tcpa(l,l,n,k)+conc(ji)/t(l,ji,k)
      do 84 l=1,mmx
      do 84 k=1,kmx
   84 tcpa(l,l,n,k)=1d0/tcpa(l,l,n,k)
      endif
   80 continue
      alloy= .not. pure
c
c---  nk1*kmx data are grouped into 'nvec' vectors, whose length are
c---  'length' except the last one whose actual length is 'lastln'.
c     In this KKR-CPA version length must be at least multiple of
c     kmx. "length=1" works quite well.
      lastln=mod(nk1*kmx-1,length)+1
      nvec=(nk1*kmx-lastln)/length+1
      do 34 k=1,kmx
      convrg(k)=.false.
   34 gfree(k)=cunit*sqrt((1d0+relfct*e(k))*e(k))
      do 36 j=1,(2*mxl-1)**2
      l=sqrt(dble(j)-5d-1)
      do 36 k=1,kmx
   36 ess(k,j)=unit*(unit/sqrt(e(k)))**l
c
      cnv=.false.
      dmpc=dmpc0
      do 40 itrcpa=1,itrmx
      if(mod(itrcpa,20) .eq. 0)dmpc=5d-1*dmpc
c---  CPA iterations will be excecuted unless cnv turns true.
      if(.not. cnv) then
      itrfin=itrcpa
c---  First, calculate LU decomposition of coherent k-matrix,
c     which will be used in the calculation of the phase.
      do 41 k=1,kmx
      if(.not. convrg(k)) then
      call equarc(tcpa(1,1,1,k),ck(1,1,1,k),mmxl**2*natm)
      do 43 i=1,natm
      mmx=(lmxtyp(itype(i))+1)**2
   43 call cprinv(ck(1,1,i,k),mmxl,mmx,1)
      call clrarc(ff(1,1,1,k),mmxl**2*natm)
      detl(k)=(0d0,0d0)
      endif
   41 continue
c---  loop over vectors composed of k-mesh and e-mesh.
      sw=0d0
      snor=0d0
!$omp parallel default(shared) private(nrun,ibsadr,gm,gi,det)
      allocate(gi(msiz**2*length),gm(msiz**2*length),det(length))
!$omp do ordered
c!$omp do
      do 42 ivec=1,nvec
      if(ivec .eq. nvec) then
      nrun=lastln
      else
      nrun=length
      endif
      ibsadr=length*(ivec-1)
c---  construct kkr matrix
      call kkrsed(gm,gi,e,tcpa,vc,unit,np,tch,pexf,prr,hh,tc,gfree
     &           ,clks,ng,ncub,iblk,itblk,its,natm,last,mxl,nd,nk
     &           ,kmx,nrun,ibsadr,ess,korder,convrg,ls,lmxblk,iatmp
     &           ,lmxtyp,itype,msiz,uu)
c     call dspgm(gm,mmxl*natm,nrun)
c---  take the inverse of the kkr matrix
      call cinvrx(gm,gi,msiz,nrun,nrun,ibsadr,kmx,convrg)
c---  phase of T-matrix is calcualted.
      call getdtb(e,det,ck,unit,gm,prr,np,itype,ibsadr,kmx,nk,mxl
     &           ,msiz,nrun,natm,korder,convrg,lmxtyp,iatmp)
c---  take k-summation
!$omp ordered
c!$omp critical
      call bzmsmb(gi,det,ff,detl,wtkp,nk,mxl,natm,kmx,ibsadr
     &     ,nrun,sw,snor,korder,convrg,lmxtyp,itype,iatmp,msiz)
c!$omp end critical
!$omp end ordered
   42 continue
      deallocate(gi,gm,det)
!$omp end parallel
      snor=snor/dble(kmx)
c---  rotate ff and sum them up.
c     The normalization ff -> ff/snor will be also taken.
      call bzmrot(ff,natm,irotat,urotat,isymop,kmx,convrg
     &           ,mxl,snor,lmxtyp,itype,ls)
c---  coherent k-matrix is iterated by a single step following the
c     Nesbet's algorithm or ATA algorithm, which normaly are stable.
c     --- Nesbet's iteration algorithm used,
      if(nesitr) then
      call nesbet(t,tcpa,ff,cj,conc,dmpc,kmx,ncmp,ntyp,ncmpx,natm
     &           ,itype,convrg,mxl,alloy,lmxtyp)
      else
c     --- ATA(average t-matrix approximation) iteration algorithm used.
      call cpaitr(t,tcpa,ff,cj,conc,dmpc,kmx,ncmp,ntyp,ncmpx,natm
     &           ,itype,convrg,mxl,alloy,lmxtyp)
      endif
      call bzmrot(tcpa,natm,irotat,urotat,isymop,kmx,convrg
     &           ,mxl,1d0,lmxtyp,itype,ls)
c---  check if coherent k-matrix converges.
      cnv=.true.
      do 48 k=1,kmx
   48 cnv=cnv .and. convrg(k)
      endif
   40 continue
c---  check the continuity in the phase of "str" and construct
c     the local green's function for each site and component.
c     do 61 k=1,kmx
c  61 write(*,'(i3,10f12.7)')k,(str(l,4,k),l=5,9)
c  61 write(*,'(i3,10f12.7)')k,t(9,4,k)
      call sbrnch(str,mmxl*ncmpx,kmx)
c     do 62 k=1,kmx
c  62 write(*,'(i3,10f12.7)')k,(str(l,4,k),l=5,9)
      call green(g,pf,phase,str,ck,ff,cj,t,mxl,ntyp,ncmp,iatm,kmx
     &          ,natm,ncmpx,lmxtyp)
c     call rotave(g,urotat,isymop,kmx,mxl,lmxtyp,ncmpx)
      call sbrnch(phase,mmxl*ncmpx,kmx)
      do 64 k=1,kmx
   64 detl(k)=detl(k)/snor
c  64 detl(k)=(0d0,0d0)
      ji=0
      do 66 i=1,ntyp
      mmx=(lmxtyp(i)+1)**2
      do 66 j=1,ncmp(i)
      ji=ji+1
      wt=conc(ji)*dble(iwtyp(i))
      do 66 l=1,mmx
c     do 66 l=5,9
      do 66 k=1,kmx
   66 detl(k)=detl(k)-wt*phase(l,ji,k)
c  66 detl(k)=detl(k)-wt*str(l,ji,k)
c  66 detl(k)=detl(k)
c     call sbrnch(detl,1,kmx)
c
      do 70 k=2,kmx
   70 detl(k)=(detl(k)-detl(1))/pi
      detl(1)=(0d0,0d0)
c     ------ print out detl -------
      if(iprint(3) .ge. 1) then
      write(*,1300)(k,detl(k),k=1,kmx)
 1300 format(/3x,'detl(k)'/(1x,i3,3x,2f12.7))
      endif
c
c     --- calculation of the Bloch spectral functions
      if(ids .eq. 4) then
      do 130 k=1,kmx
  130 convr2(k)=.false.
      inivec=nvec
      if(lastln .eq. length) inivec=nvec+1
      lastln=mod(nk*kmx-1,length)+1
      nvec=(nk*kmx-lastln)/length+1
!$omp parallel default(shared) private(nrun,ibsadr,gm,gi,det)
      allocate(gi(msiz**2*length),gm(msiz**2*length),det(length))
!$omp do
      do 100 ivec=inivec,nvec
      if(ivec .eq. nvec) then
      nrun=lastln
      else
      nrun=length
      endif
      ibsadr=length*(ivec-1)
c---  construct kkr matrix
      call kkrsed(gm,gi,e,tcpa,vc,unit,np,tch,pexf,prr,hh,tc,gfree
     &           ,clks,ng,ncub,iblk,itblk,its,natm,last,mxl,nd,nk
     &           ,kmx,nrun,ibsadr,ess,korder,convr2,ls,lmxblk,iatmp
     &           ,lmxtyp,itype,msiz,uu)
c---  take the inverse of the kkr matrix
      call cinvrx(gm,gi,msiz,nrun,nrun,ibsadr,kmx,convr2)
c---  phase of T-matrix is calcualted.
      call getdtb(e,det,ck,unit,gm,prr,np,itype,ibsadr,kmx,nk,mxl
     &           ,msiz,nrun,natm,korder,convr2,lmxtyp,iatmp)
      do 100 k=1,nrun
      iadr=ibsadr+k
      kadr=(iadr-1)/kmx-nk1+1
      nadr=mod(iadr-1,kmx)+1
  100 if(kadr .ge. 1) spctrl(nadr,kadr)=-det(k)
      deallocate(gi,gm,det)
!$omp end parallel
      call clrarc(gfree,kmx)
      ji=0
      do 110 i=1,ntyp
      mmx=(lmxtyp(i)+1)**2
      do 110 j=1,ncmp(i)
      ji=ji+1
      wt=conc(ji)*dble(iwtyp(i))
      do 110 l=1,mmx
      do 110 k=1,kmx
c      --- gfree is used only as a work space, having nothing to do
c          with the free Green function.
  110 gfree(k)=gfree(k)+wt*phase(l,ji,k)
      do 120 k=1,kmx
      do 120 kp=1,nk3
  120 spctrl(k,kp)=(spctrl(k,kp)-gfree(k))/pi
      deallocate(convr2)
      endif
c     stop
      end
