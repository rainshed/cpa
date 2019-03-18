      subroutine fxspin(total,cnutr,anclr,anc,ef,dosef,def,sftef,sm
     &        ,ew,ez,itype,mmxl,mxlcmp,ng,natm,ntyp,ncmp,ncmpx,conc
     &        ,fspin)
c----------------------------------------------------------------------
c     -----------------------------------------
c     --- KKR-CPA fixed spin moment version ---
c     -----------------------------------------
c     Check charge neutralityi and the spin moment.
c     If they differ from the expected ones, shift the Fermi levels
c     for up and down spins separately.
c     Also the corrections for sm are made.
c     coded by H.Akai, 1987, Osaka
c     KKR-CPA implemented by H.Akai, 19 Sep. 1996, Osaka
c     Fixed spin moment procedure implemented by H. Akai 5. Aug. 99
c     modified by H. Akai, 25 Aug. 1999, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension anclr(ncmpx),anc(ncmpx),total(2),tr(30)
     &         ,dosef(mmxl,ncmpx,2),ef(2),def(2)
     &         ,sm(ng,mmxl,ncmpx,2),conc(ncmpx),sftef(2)
      integer itype(natm),ncmp(ntyp),mxlcmp(ncmpx)
      data slimit/5d-2/, amrgin/9d-1/
      cnutr=total(1)+total(2)
      do 10 n=1,natm
      i=itype(n)
      do 10 j=1,ncmp(i)
      call jip(i,j,ji)
   10 cnutr=cnutr+conc(ji)*(dble(int(anc(ji)+amrgin))-anclr(ji))
      tspin=total(1)-total(2)-fspin
      sftef(1)=5d-1*(cnutr+tspin)
      sftef(2)=5d-1*(cnutr-tspin)
      do 20 is=1,2
      sftef(is)=-sftef(is)/def(is)
      sftef(is)=min(slimit,max(-slimit,sftef(is)))
   20 total(is)=total(is)+sftef(is)*def(is)
      do 30 is=1,2
      call gntcs(ef(is),ew,ez,tr,ng)
      do 30 i=1,ncmpx
      do 30 l=1,mxlcmp(i)**2
      d=dosef(l,i,is)*sftef(is)
      do 30 n=1,ng
   30 sm(n,l,i,is)=sm(n,l,i,is)+d*tr(n)
      return
      end
