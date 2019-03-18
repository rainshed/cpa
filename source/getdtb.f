      subroutine getdtb(e,det,ck,unit,gm,prr,np,itype,ibsadr,kmx,nk,mxl
     &                 ,msiz,nrun,natm,korder,convrg,lmxtyp,iatmp)
c-----------------------------------------------------------------------
c     -----------------------
c     --- KKR-CPA version ---
c     -----------------------
c     Given a Green function of trianguler matrix form 'gm', this
c     returns the logarithm of the determinant of the Green function.
c     This originary is a tiny part of 'conkkr'.
c     Coded by H.Akai, 11 Feb. 1996, Osaka
c
c     Only difference from the usual KKR version is that now the
c     k-matrix has mxl**2 elements instead of mxl.
c     KKR CPA implemented by H.Akai, 15 Sep 1996, Osaka
c     Modified by H. Akai, 25 Aug. 1999, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 e(kmx),det(nrun),ck(mxl**2,mxl**2,natm,kmx)
     &          ,gm(msiz,msiz,nrun),add
      real*8 prr(np,nk)
      integer itype(natm),korder(nk),iatmp(natm),lmxtyp(*)
      logical convrg(kmx)
      call sbtime(3,0)
      eunit=unit**2
c     ---back scattering contribution
      do 10 k=1,nrun
      iadr=ibsadr+k
      kadr=korder((iadr-1)/kmx+1)
      nadr=mod(iadr-1,kmx)+1
      if(.not. convrg(nadr)) then
      det(k)=(0d0,0d0)
      do 20 n=1,natm
      iap=iatmp(n)-1
      mmx=(lmxtyp(itype(n))+1)**2
      do 20 j=1,mmx
      jj=j+iap
      add=gm(jj,jj,k)/ck(j,j,n,nadr)
   20 det(k)=det(k)+log(add)
c
c     ---free branch contribution
      do 30 i=1,np
      add=prr(i,kadr)-e(nadr)/eunit
   30 det(k)=det(k)+log(add)
      endif
   10 continue
c---   The following codes may be used in order to speed up
c      the phase calculation by avoiding log.
c      add=gm(jj,jj,k)/ck(j,j,n,nadr)
c      sd=sign(1d0,dimag(det(ik)))
c      det(k)=det(k)*add
c      br(k)=br(k)+sa*2.5d-1*(sd+sa)*(sd-sign(1d0,dimag(det(k))))
c
c      add=prr(i,kadr)-e(nadr)/eunit
c      sd=sign(1d0,dimag(det(k)))
c      sa=sign(1d0,dimag(add))
c      det(k)=det(k)*add
c      br(k)=br(k)+sa*2.5d-1*(sd+sa)*(sd-sign(1d0,dimag(det(k))))
      call sbtime(3,1)
      end
