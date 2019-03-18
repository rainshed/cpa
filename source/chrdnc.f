      subroutine chrdnc(ro,rorg,tof,f,mxj,mxlcmp,sm,tm,ng,ntyp,xr
     &                 ,meshr,rstr,total,itype,natm,ncmp,ncmpx,conc)
c----------------------------------------------------------------------
c     ------------------------------------
c     --- KKR-CPA + spin-orbit version ---
c     ------------------------------------
c     Construct charge density from sm data combined with radial
c     wave function data.
c     coded by H.Akai, 1983, Juelich
c     KKR-CPA implemented by H.Akai, 19 Sep. 1996, Osaka
c     adapted for spin-orbit version, 18 April, 1997, Osaka
c     modified by H. Akai, 25 Aug. 1999, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension tm(ng,ng),ro(meshr,ncmpx),rorg(20,ncmpx)
     &         ,xr(meshr,ncmpx),rstr(meshr,mxj,ng,ncmpx)
     &         ,tof(mxj,ncmpx),f(mxj,ncmpx,ng)
     &         ,sm(ng,mxj,ncmpx),conc(ncmpx)
      integer itype(natm),ncmp(ntyp),mxlcmp(ncmpx)
      call reduce(sm,f,tm,ng,mxj*ncmpx)
      call clrarr(tof,ncmpx*mxj)
      do 10 n=1,ng
      do 10 i=1,ncmpx
      do 10 l=1,mxlcmp(i)**2
      tof(l,i)=tof(l,i)+f(l,i,n)
c     write(*,'(3i3,1pe15.7)')n,i,l,f(l,i,n)
      do 10 k=1,meshr-1
   10 ro(k,i)=ro(k,i)+rstr(k,l,n,i)**2*f(l,i,n)
      do 20 i=1,ncmpx
      do 30 k=1,meshr-1
   30 ro(k,i)=ro(k,i)/xr(k,i)**2
   20 call extorg(rorg(20,i),ro(1,i),xr(1,i))
      roi=total
      omi=0d0
      do 80 n=1,natm
      i=itype(n)
      do 80 j=1,ncmp(i)
      call jip(i,j,ji)
      omi=omi+conc(ji)*(xr(meshr,ji)**3-xr(meshr-1,ji)**3)
      roi=roi+conc(ji)*ro(meshr,ji)
      do 80 l=1,mxlcmp(ji)**2
   80 roi=roi-conc(ji)*tof(l,ji)
      do 90 i=1,ncmpx
   90 ro(meshr,i)=3d0*roi/omi
c     write(*,'(1x,1p,16e13.5)')(f(l,1,1),l=1,16)
      end
