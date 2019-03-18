      subroutine neutrl(total,cnutr,anclr,anc,ef,dosef,def,sftef,sm
     &          ,ew,ez,itype,mmxl,mxlcmp,ng,natm,ntyp,ncmp,ncmpx,conc)
c----------------------------------------------------------------------
c     -----------------------
c     --- KKR-CPA version ---
c     -----------------------
c     Check charge neutrality. If not, shift the fermi level
c     and make correction for sm to follow it up.
c     coded by H.Akai, 1987, Osaka
c     KKR-CPA implemented by H.Akai, 19 Sep. 1996, Osaka
c     modified by H. Akai, 25 Aug. 1999, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension anclr(ncmpx),anc(ncmpx),total(2),tr(30)
     &         ,dosef(mmxl,ncmpx,2),ef(2),def(2)
     &         ,sm(ng,mmxl,ncmpx,2),conc(ncmpx)
      integer itype(natm),ncmp(ntyp),mxlcmp(ncmpx)
      data slimit/5d-2/ , amrgin/9d-1/
      cnutr=total(1)+total(2)
      do 10 n=1,natm
      i=itype(n)
      do 10 j=1,ncmp(i)
      call jip(i,j,ji)
   10 cnutr=cnutr+conc(ji)*(dble(int(anc(ji)+amrgin))-anclr(ji))
      sftef=-cnutr/(def(1)+def(2))
      sftef=min(slimit,max(-slimit,sftef))
      total(1)=total(1)+sftef*def(1)
      total(2)=total(2)+sftef*def(2)
      do 20 is=1,2
      call gntcs(ef(is),ew,ez,tr,ng)
      do 20 i=1,ncmpx
      do 20 l=1,mxlcmp(i)**2
      d=dosef(l,i,is)*sftef
      do 20 n=1,ng
   20 sm(n,l,i,is)=sm(n,l,i,is)+d*tr(n)
      end
