      subroutine potenv2(sdftyp,itype,ntyp,natm,anclr,ro,v,q,exspl
     &                 ,dr,xr,meshr,a,w,amdlng,u,ncmp,ncmpx,conc
     &                 ,efield,atmicp,iatm)
c----------------------------------------------------------------------
c     -----------------------
c     --- KKR-CPA version ---
c     -----------------------
c     +---------------------------------------------------+
c        sdf parametrization by vbh     ... vbh
C                            by mjw     ... mjw
c                            by vwn     ... vwn
c     Langreth, Mehl non-local + mjw    ... lmmjw
c                   non-local + vbh     ... lmvbh  missing
c     Perdew, Yue    non-local + mjw    ... pymjw
c                   non-local + vbh     ... pyvbh  missing
c                   non-local + vwn     ... pyvwn
C     Perdew, Wang   GGA91              ... gga91
c     Engel, Vosko   GGA exchange-only  ... ev
c     +---------------------------------------------------+
c     Generate potentials from charge/spin densities.
c
c     coded by H.Akai, 1984, Juelich
c     modified to include gradient terms by H.Ebert, 1991, Erlangen
c     modified to fit complex lattices by H.Akai, 1992, Osaka
c     KKR-CPA implemented by H.Akai, 13 Sep, 1996, Osaka
c     Further implementation of GGA M. Batocretti, Muenchen
c     Adapted to cpa97 and cpa98 coded H. Akai, August, 1998, Duisburg
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      real*8 ro(meshr,ncmpx,2),v(meshr,ncmpx,2),dr(meshr,ncmpx)
     &      ,xr(meshr,ncmpx),anclr(ncmpx),w(meshr,ncmpx),q(ncmpx)
     &      ,amdlng(natm,natm),conc(ncmpx),atmicp(3,natm)
      integer itype(natm),ncmp(ntyp),iatm(ntyp)
      character*6 sdftyp
c --- w(k,i,1): hartree + xc energy at k
c --- w(k,i,2): u - 3(exc-vexc)*roh, used for 3pv calculation.<--missing
      pi=4d0*atan(1d0)
      do 10 i=1,ncmpx
      call poisna(ro(1,i,1),v(1,i,1),anclr(i),dr(1,i),xr(1,i),meshr)
   10 q(i)=v(meshr-1,i,1)*xr(meshr-1,i)/2d0
c  10 write(*,*)'i,q(i)=',i,q(i)
c     --- qint is the interstitial charge density.
c     --- vint=interstitial volume/4pi
      qint=0d0
      vint=0d0
      do 20 ii=1,natm
      i=itype(ii)
      do 20 j=1,ncmp(i)
      call jip(i,j,ji)
      qint=qint-q(ji)*conc(ji)
   20 vint=vint+(xr(meshr,ji)**3-xr(meshr-1,ji)**3)*conc(ji)
      if(vint .gt. 1d-10) then
      vint=vint/3d0
      qint=qint/vint
      else
      vint=0d0
      qint=0d0
      endif
c     --- charge density at the interstitial region is set 'qint'
c         such that it satisfies the cahrge neutrality condition
c         in the unit cell irrespective of the actual charge
c         neutrality obtained with trial fermi level. the spin
c         density on the other hand is kept unchanged.
      ji=0
      do 30 i=1,ntyp
      do 30 j=1,ncmp(i)
      ji=ji+1
      do 40 k=1,natm
      kk=k
      if(i .eq. itype(k)) go to 50
   40 continue
      call errtrp(1,'potenv','type not found')
   50 ro(meshr,ji,1)=qint
      v(meshr,ji,1)=0d0
      do 60 ia=1,natm
      l=itype(ia)
      c=amdlng(kk,ia)/a
      do 60 jj=1,ncmp(l)
      call jip(l,jj,jjl)
   60 v(meshr,ji,1)=v(meshr,ji,1)+c*q(jjl)*conc(jjl)
      do 70 k=1,meshr-1
   70 w(k,ji)=5d-1*v(k,ji,1)-anclr(ji)/xr(k,ji)
      w(meshr,ji)=0d0
   30 call equarr(v(1,ji,1),v(1,ji,2),meshr)
      u=0d0
      do 80 ia=1,natm
      i=itype(ia)
      do 80 j=1,ncmp(i)
      call jip(i,j,ji)
   80 u=u-conc(ji)*q(ji)*v(meshr,ji,1)/2d0
c
c     --- now choose a type of xc potential.
      if(sdftyp .eq. 'vbh') then
      call excvbh(ro,v,w,ncmpx*meshr,ncmpx*meshr)
      elseif(sdftyp .eq. 'mjw') then
      call excmjw(ro,v,w,ncmpx*meshr,ncmpx*meshr)
      elseif(sdftyp .eq. 'vwn') then
      call excvwn(ro,v,w,ncmpx*meshr,ncmpx*meshr)
      elseif(sdftyp .eq. 'lmmjw') then
      call exclmm(ro,v,w,dr,xr,meshr,ncmpx,meshr)
      elseif(sdftyp .eq. 'pymjw') then
      call excpym(ro,v,w,dr,xr,meshr,ncmpx,meshr)
      elseif(sdftyp .eq. 'pyvwn') then
      call excpyv(ro,v,w,dr,xr,meshr,ncmpx,meshr)
      elseif(sdftyp .eq. 'gga91') then
      call excg91(ro,v,w,dr,xr,meshr,ncmpx,meshr)
      elseif(sdftyp .eq. 'ev') then
      call excev(ro,v,w,dr,xr,meshr,ncmpx,meshr)
c     elseif(sdftyp .eq. 'pbe') then
c     call excpbe(ro,v,w,dr,xr,meshr,ncmpx,meshr)
      else
      call errtrp(1,'potsic','sdftyp '//sdftyp//' not serviced')
      endif
      exspl=v(meshr,1,2)-v(meshr,1,1)
      do 90 is=1,2
      do 90 i=1,ntyp
c     --- electric field along z-direction is applied.
c     anode and cathod are placed at -z and z direction respectively.
      bias=efield*atmicp(3,iatm(i))
      do 90 j=1,ncmp(i)
      call jip(i,j,ji)
      do 90 k=1,meshr
   90 v(k,ji,is)=v(k,ji,is)-v(meshr,ji,is)+bias
      u=u+w(meshr,1)*vint
      do 100 i=1,ncmpx
  100 w(meshr,i)=fintgr(w(1,i),dr(1,i),xr(1,i),meshr)
      do 110 ia=1,natm
      i=itype(ia)
      do 110 j=1,ncmp(i)
      call jip(i,j,ji)
  110 u=u+conc(ji)*w(meshr,ji)
      end
