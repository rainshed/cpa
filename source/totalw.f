      subroutine totalw(ro,v,w,corlvl,bnd2,esic,meshr,xr,dr
     &                 ,natm,ntyp,itype,u,te,config,ncmp,ncmpx,conc)
c----------------------------------------------------------------------
c     -----------------------
c     --- KKR-CPA version ---
c     -----------------------
c     Calculate total energy from charge density and input potential.
c     Hartree+xc part must have been calculated by 'potenv' or any
c     equivalnet routine and the result is given as 'u'.
c     coded by H.Akai, 1987, Osaka
c     KKR-CPA implemented by H.Akai, 21 Sep. 1996, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (ncmpmx=40)
      real*8 ro(meshr,ncmpx,2),v(meshr,ncmpx,2),w(meshr)
     &      ,corlvl(18,ncmpx,2),bnd2(2)
     &      ,xr(meshr,ncmpx),dr(meshr,ncmpx)
     &      ,config(18,ncmpx,2),esic(ncmpx,2)
     &      ,conc(ncmpx)
      real*8,allocatable::t(:)
      integer itype(natm),ncmp(ntyp)
      allocate(t(ncmpx))
      te=u+bnd2(1)+bnd2(2)
      do 10 i=1,ncmpx
      do 20 k=1,meshr-1
   20 w(k)=((ro(k,i,1)+ro(k,i,2))*5d-1*(v(k,i,1)-v(meshr,i,1))
     &     +(ro(k,i,1)-ro(k,i,2))*5d-1*(v(k,i,2)-v(meshr,i,2))
     &         )*xr(k,i)**2*dr(k,i)
      t(i)=0d0
      do 30 k=1,meshr-3,2
   30 t(i)=t(i)+w(k)+4d0*w(k+1)+w(k+2)
      t(i)=-t(i)/3d0+esic(i,1)+esic(i,2)
      do 10 is=1,2
      do 10 n=1,18
   10 t(i)=t(i)+corlvl(n,i,is)*max(0d0,config(n,i,is))
      do 40 n=1,natm
      i=itype(n)
      do 40 j=1,ncmp(i)
      call jip(i,j,ji)
   40 te=te+conc(ji)*t(ji)
      end
