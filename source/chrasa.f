      subroutine chrasa(ro,rorg,tof,f,mxj,mxlcmp,sm,tm,ng,ntyp,xr,meshr
     &                 ,rstr,total,itype,natm,ncmp,ncmpx,conc,lmxtyp)
c----------------------------------------------------------------------
c     Construct charge density from sm data combined with radial
c     wave function data.
c     coded by H.Akai, 1983, Juelich
c     This version is used for KKR-ASA. The charge density inside
c     the sphere is renormalized so as to give the correct number
c     of electrons in the system.
c     coded by H.Akai, 1995, Duisburg
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension tm(ng,ng),ro(meshr,ncmpx),rorg(20,ncmpx)
     &         ,xr(meshr,ncmpx),rstr(meshr,mxj,ng,ncmpx)
     &         ,tof(mxj,ncmpx),f(mxj,ncmpx,ng)
     &         ,sm(ng,mxj,ncmpx),conc(ncmpx)
      integer itype(natm),ncmp(ntyp),mxlcmp(ncmpx),lmxtyp(ntyp)
      call reduce(sm,f,tm,ng,mxj*ncmpx)
      call clrarr(tof,ncmpx*mxj)
      totin=0d0
      do 50 i=1,natm
      j=itype(i)
      mmx=(lmxtyp(j)+1)**2
      do 50 in=1,ncmp(j)
      call jip(j,in,jn)
      do 50 n=1,ng
      do 50 l=1,mmx
   50 totin=totin+conc(jn)*f(l,jn,n)
c     write(*,*)'totin=',totin,'  total=',total
      renorm=total/totin
      do 10 n=1,ng
      do 10 i=1,ncmpx
      do 10 l=1,mxlcmp(i)**2
      f(l,i,n)=renorm*f(l,i,n)
      tof(l,i)=tof(l,i)+f(l,i,n)
      do 10 k=1,meshr-1
   10 ro(k,i)=ro(k,i)+rstr(k,l,n,i)**2*f(l,i,n)
      do 20 i=1,ncmpx
      do 30 k=1,meshr-1
   30 ro(k,i)=ro(k,i)/xr(k,i)**2
   20 call extorg(rorg(20,i),ro(1,i),xr(1,i))
      do 40 i=1,ncmpx
   40 ro(meshr,i)=0d0
      end
