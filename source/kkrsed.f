      subroutine kkrsed(gm,d,e,t,vc,unit,np,tch,pexf,prr,hh,tc,gfree
     &                 ,clks,ng,ncub,iblk,itblk,its,natm,last,mxl
     &                 ,nd,nk,kmx,nrun,ibsadr,ess,korder,convrg,ls
     &                 ,lmxblk,iatmp,lmxtyp,itype,msiz,uu)
c-----------------------------------------------------------------------
c     In this version the phase factor i^(l1-l2) is included
c     in the Green's function.
c     Modified by H. Akai, 30 July 2005, Osaka
c
c     ------------------------------------
c     --- KKR-CPA + spin-orbit version ---
c         + variable l_max
c     ------------------------------------
c     Calculate structure constants from their tchebycheff expansion
c     coefficient and construct kkr matrix.
c     Originally coded as 'dless' by M.Akai, 1980, Osaka,
c     revised as 'prset' by M.Akai and H.Akai, 1981, Osaka,
c     revised as 'kkrset' by H.Akai, 1986, Juelich,
c     latest revision by H.Akai, Feb. 1996, Osaka.
c
c     Main modification is now k-point is chosen according to
c     a randomized order stored in 'korder'.
c     Another modification is now t-matrix has full diagonal
c     elements, t(1),t(2),...t(mxl**2).
c     KKR CPA implemented by H.Akai, 14 Sep 1996, Osaka.
c     Adapted to spin-orbit version, 25 April, 1997, Osaka.
c     Modified so as adapted to the variable lmax cases
c     by H. Akai, 25 Aug. 1999, Osaka.
c-----------------------------------------------------------------------
c     --- t is the inverse of the negative t-matrix ---
      implicit real*8 (a-h,o-z)
      complex*16 gm(msiz,msiz,nrun)
     &          ,t(mxl**2,mxl**2,natm,kmx),d((2*mxl-1)**2,nd)
     &          ,tc(ng,kmx),tch((2*mxl-1)**2,ng,nk,nd)
     &          ,pexf(np,nk,nd),e(kmx),gfree(kmx)
     &          ,ess(kmx,(2*mxl-1)**2),eu,phf,phfcl
     &          ,uu(2*mxl-1,2*mxl-1,2)
      real*8 prr(np,nk),hh(np,(2*mxl-1)**2,nk),clks(last)
      integer ncub(last),iblk(natm,natm),itblk(5,its),itype(natm)
     &       ,korder(nk),lmxblk(nd),iatmp(natm),lmxtyp(*)
      logical convrg(kmx)
      call sbtime(1,0)
      pi=4d0*atan(1d0)
      q1=2d0/pi/vc
      mmxl=mxl**2
      mmxj=(2*mxl-1)**2
c     basadr=dble(ibsadr)
c     enblk=dble(kmx)
c     ---relativistic treatment here is not very correct: only
c     ---the correction associated to the single site scattering term
c     ---is considered.
      eunit=unit**2
      do 140 k=1,nrun
      iadr=ibsadr+k
      kadr=korder((iadr-1)/kmx+1)
      nadr=mod(iadr-1,kmx)+1
      if(.not. convrg(nadr)) then
      eu=e(nadr)/eunit
      do 40 i=1,nd
      mmx=(lmxblk(i)+1)**2
      do 40 mr=1,mmx
   40 d(mr,i)=(0d0,0d0)
      do 50 n=1,np
      do 50 i=1,nd
      mmx=(lmxblk(i)+1)**2
      do 50 mr=1,mmx
   50 d(mr,i)=d(mr,i)
     &   +pexf(n,kadr,i)*hh(n,mr,kadr)/(eu-prr(n,kadr))
      do 60 i=1,nd
      mmx=(lmxblk(i)+1)**2
      do 60 mr=1,mmx
   60 d(mr,i)=q1*d(mr,i)
c
      do 70 n=1,ng
      do 70 i=1,nd
      mmx=(lmxblk(i)+1)**2
      do 70 mr=1,mmx
   70 d(mr,i)=d(mr,i)+tch(mr,n,kadr,i)*tc(n,nadr)
      do 80 i=1,nd
      mmx=(lmxblk(i)+1)**2
      do 80 mr=1,mmx
   80 d(mr,i)=d(mr,i)*ess(nadr,mr)
c
      do 100 j=1,msiz
      do 100 i=1,msiz
  100 gm(i,j,k)=(0d0,0d0)
      ji=1
      do 90 mc=1,mmxl
      l2=sqrt(dble(mc)-5d-1)
      do 90 mr=1,mmxl
      l1=sqrt(dble(mr)-5d-1)
c     --- If a phase factor should be attached,activate the following.
      phf=(0d0,1d0)**(l1-l2)
c     --- Otherwise, use the following.
c     phf=(1d0,0d0)**(l1-l2)
c     --- In this case, the phase factor should be taken into account
c         later. This affect the procedure taking the inversion
c         k -> -k in the subroutine bzmrot. Also the
c         phase factor should be attached 'phf' calculated in the
c         subroutine cpshft.
c
c     phf=1d0
      do 110 j=ji,last
      jj=j
      m=ncub(j)
      if(m .eq. 0) go to 90
      phfcl=clks(j)*phf
      do 110 i=1,its
      if(mr .le. (itblk(2,i)+1)**2) then
      if(mc .le. (itblk(3,i)+1)**2) then
      iap=iatmp(itblk(4,i))+mr-1
      ibp=iatmp(itblk(5,i))+mc-1
      gm(iap,ibp,k)=gm(iap,ibp,k)+d(m,itblk(1,i))*phfcl
c     gm(iap,ibp,k)=gm(iap,ibp,k)+d(m,itblk(1,i))*clks(j)
      endif
      endif
  110 continue
   90 ji=jj+1
      if(ls .eq. 1) then
      do 180 i=1,its
      mxla=itblk(2,i)+1
      mxlb=itblk(3,i)+1
      iap=iatmp(itblk(4,i))
      ibp=iatmp(itblk(5,i))
  180 call gmtrns(gm(1,1,k),1,mxla,mxlb,msiz,iap,ibp,'r2s',uu,mxl)
      endif
      do 120 ib=1,natm
      do 120 ia=1,natm
      i=iblk(ia,ib)
      ka=itblk(4,i)
      kb=itblk(5,i)
c     if(nadr .eq. kmx) then
c     write(*,'('' ia,ib'',2i3,'' i,ka,kb'',3i3)')ia,ib,i,ka,kb
c     endif
      if(ia .ne. ka .or. ib .ne. kb) then
      iap=iatmp(ia)-1
      ibp=iatmp(ib)-1
      kap=iatmp(ka)-1
      kbp=iatmp(kb)-1
      mmxa=(lmxtyp(itype(ia))+1)**2
      mmxb=(lmxtyp(itype(ib))+1)**2
c     if(nadr .eq. kmx) then
c     write(*,'('' iap,ibp,kap,kbp,mmxa,mmxb'',4i3,2i4)')
c    & iap,ibp,kap,kbp,mmxa,mmxb
c     endif
      do 130 mc=1,mmxb
      do 130 mr=1,mmxa
  130 gm(mr+iap,mc+ibp,k)=gm(mr+kap,mc+kbp,k)
      endif
  120 continue
c     if(nadr .eq. kmx) then
c     do 160 ib=1,natm
c     do 160 ia=1,natm
c     i=iblk(ia,ib)
c     iap=iatmp(ia)-1
c     ibp=iatmp(ib)-1
c     write(*,'(1x,''block='',i3,'' atoma, atomb'',2i3)')i,ia,ib
c 160 write(*,'(1x,1p,6e13.6)')((gm(mr+iap,mc+ibp,k),mr=1,1),mc=1,1)
c     endif
      do 170 ia=1,natm
      iap=iatmp(ia)-1
      mmx=(lmxtyp(itype(ia))+1)**2
      do 170 mc=1,mmx
      l2=sqrt(dble(mc)-5d-1)
      gm(mc+iap,mc+iap,k)=gm(mc+iap,mc+iap,k)+gfree(nadr)
      do 170 mr=1,mmx
      l1=sqrt(dble(mr)-5d-1)
c 170 gm(mr+iap,mc+iap,k)=gm(mr+iap,mc+iap,k)*(0d0,1d0)**(l1-l2)
c    &                   +t(mr,mc,ia,nadr)
  170 gm(mr+iap,mc+iap,k)=gm(mr+iap,mc+iap,k)+t(mr,mc,ia,nadr)
c     do 170 ia=1,natm
c     do 170 ib=1,natm
c     if(ia .eq. ib) then
c     do 190 m=1,mmxl
c 190 gm(m,ia,m,ia,k)=gm(m,ia,m,ia,k)+gfree(nadr)
c     endif
c     do 200 mr=1,mmxl
c     do 210 mc=1,mmxl
c     d(mc,1)=set
c     do 210 m=1,mmxl
c 210 d(mc,1)=d(mc,1)+gm(mr,ia,m,ib,k)*t(m,mc,ib,nadr)
c     if(ia .eq. ib) d(mr,1)=d(mr,1)+1d0
c     do 200 mc=1,mmxl
c 200 gm(mr,ia,mc,ib,k)=d(mc,1)
c 170 continue
      endif
c 140 continue
c     write(*,*)'iadr=',iadr,'  k',k
  140 continue
c     stop
c     endif
c 150 continue
      call sbtime(1,1)
      return
      end
