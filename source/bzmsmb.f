      subroutine bzmsmb(gi,det,ff,detl,wtkp,nk,mxl,natm,kmx,ibsadr
     &           ,nrun,sw,snor,korder,convrg,lmxtyp,itype,iatmp,msiz)
c---------------------------------------------------------------------
c     -----------------------
c     --- KKR-CPA version ---
c     -----------------------
c     Note: valid only when length=me, where
c     'length' and 'me' is defined and used in the calling routine.
c     Integrate the Green function on BZ mesh.
c     coded by H.Akai, 1989, Osaka
c     latest revision, 11 Feb. 1996, Osaka
c     KKR CPA implemented by H.Akai, 14 Sep 1996, Osaka
c     A small modification (the oeder of the index of ff is
c     changed), 13 August, 1997, Duisburg
c     Modifeid by H. Akai, 25 Aug. 1999, Osaka
c     Now adapted to arbitrary "length" used in the calling routine.
c     Modified by H. Akai, 5 Jan, 2016, Tokyo
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 gi(msiz,msiz,nrun)
     &          ,det(nrun),ff(mxl**2,mxl**2,natm,kmx),detl(kmx)
      real*8 wtkp(nk)
      integer korder(nk),lmxtyp(*),itype(natm),iatmp(natm)
      logical convrg(kmx)
      do 20 k=1,nrun
      iadr=ibsadr+k
      kadr=korder((iadr-1)/kmx+1)
      rdwt=wtkp(kadr)
      nadr=mod(iadr-1,kmx)+1
c     --- The weight 'rdwt' is summed up every time the loop 20 is
c         executed. Therefore the final 'snor' is kmx times more
c         than the actual 'snor'. This has to be further corrected
c         after the energy*k-point loop ends in the calling routine.
      snor=snor+rdwt
      if(.not. convrg(nadr)) then
      do 10 i=1,natm
      mmx=(lmxtyp(itype(i))+1)**2
      iap=iatmp(i)-1
      do 10 mc=1,mmx
      do 10 mr=1,mmx
c     ---add up the site-diagonal block of gi to obtain f(k,l1,l2,i).
   10 ff(mr,mc,i,nadr)=ff(mr,mc,i,nadr)+rdwt*gi(mr+iap,mc+iap,k)
      detl(nadr)=detl(nadr)-rdwt*det(k)
      endif
   20 continue
c     iadr=ibsadr+nrun
c     nadr=mod(iadr-1,kmx)+1
c     if(.not. convrg(nadr)) then
c     do 30 i=2,2
c     mmx=(lmxtyp(itype(i))+1)**2
c     iap=iatmp(i)-1
c     write(*,'(/a,i3,a,i2)')'bzmsmb k-point=',kadr,'  atom=',i
c     write(*,'((1x,1p9e20.13))')
c    &   ((dble(gi(mr+iap,mc+iap,nrun)),mc=1,9),mr=1,9)
c  30 write(*,'(/(1x,1p9e20.13))')
c    &   ((imag(gi(mr+iap,mc+iap,nrun)),mc=1,9),mr=1,9)
c     endif
      end
