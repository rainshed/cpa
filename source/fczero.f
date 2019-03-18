      subroutine fczero(cm,ew,ez,tchc,tchs,fcs,tm,elvl
     &                 ,ng,jmx,jmxcmp,xlim)
c-----------------------------------------------------------------------
c     This program searches all zeros of C and calculates the
c     contribution of the suprious poles of the pseudo Green's function.
c     The pseudo Green's function is defined as
c
c     G(E)=( -E**l/(C**2+E**(2*l+1)*S**2) * ((S/C)*E**(l+1)+i*sqrt(E))
c            + f/C ) * R(E)*R(E)
c         =(  i*E**(l+1/2)/(C+i*E**(l+1/2)*S) + f )*(1/C) * R(E)*R(E),
c
c     where l=0,1,2,..., and f(E) is a polynomial which fits 1/S(E)
c     at all E=En's satisfing C(En)=0. In the output 'cm' is the
c     Chebyshev expansion coefficients of f(E).
c     coded by H.Akai, 1983, Juelich
c     major modification by H.Akai, 1993, Osaka
c     very minor modification by H. Akai, 25 Aug. 1999, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 tchc(ng,jmx),tchs(ng,jmx),fcs(3,jmx),cm(ng,jmx)
     &         ,tm(ng,ng),elvl(ng)
      real*8,allocatable::d(:)
      complex*16 cx,cpolin
      complex*16,allocatable::roots(:),tc(:),y(:)
      data zero/1d-10/, small/1d-7/, t/1d-1/
c     data zero/1d-7/, small/1d-7/, t/1d-1/
      allocate(roots(ng),tc(ng),y(ng),d(ng+1))
      nr=0
      do 70 i=1,ng
      do 70 l=1,jmx
   70 cm(i,l)=0d0
      do 10 j=1,jmxcmp
      call chebpc(tchc(1,j),d,ng)
c     ---working on scaled energy centered on 'ew', scaled by 'ez'.
      e0=(fcs(3,j)-ew)/ez
      d(ng+1)=d(ng)
      do 40 i=ng,2,-1
   40 d(i)=d(i-1)-d(i)*e0
      d(1)=fcs(1,j)/ez-d(1)*e0
      big=zero
      do 110 i=1,ng+1
      big=max(big,abs(d(i)))
      if(abs(d(i)/big) .gt. zero) then
      nr0=i-1
      else
      d(i)=0d0
      endif
  110 continue
      do 120 itry=nr0,2,-1
      nr=itry
      call zroots(d,nr,roots,ierr)
      if(ierr .eq. 0) go to 130
  120 continue
      call errtrp(1,'fczero','root finding fails')
  130 do 100 i=2,nr
      if(abs(roots(i)-roots(i-1)) .lt. small) then
c     write(*,'(1x,i3,2f15.7,3x,2f15.7)')i,roots(i),roots(i-1)
      call errtrp(1,'fczero','same zero counted twice')
      endif
  100 continue
c     write(*,'(4(1p,2e13.6,2x))')(roots(i),i=1,nr)
      do 20 i=1,nr
      arg=(abs(roots(i))-xlim)/t
c     write(*,*)'arg=',arg
      arg=max(min(46d0,arg),-46d0)
      fd=1d0/(1d0+exp(arg))
c     ---now scale them back.
      roots(i)=ew+roots(i)*ez
c     write(*,'(1x,a,2i2,a,2f10.5)') 'j,i=',j,i,' z=',roots(i)
      call cgntcs(roots(i),ew,ez,tc,ng)
      y(i)=(0d0,0d0)
      do 30 l=1,ng
   30 y(i)=y(i)+tchs(l,j)*tc(l)
   20 y(i)=fd/(y(i)*(roots(i)-fcs(3,j))+fcs(2,j))
      do 60 k=1,ng
      cx=dcmplx(elvl(k),0d0)
      f=dble(cpolin(cx,roots,y,nr,err))
      do 60 l=1,ng
   60 cm(l,j)=cm(l,j)+f*tm(l,k)
c     write(*,'(4(1p,2e13.6,2x))')(roots(i),i=1,nr)
c      --- check the accuracy of the fitting.
c     do 80 i=1,nr
c     call cgntcs(roots(i),ew,ez,tc,ng)
c     cx=(0d0,0d0)
c     do 90 l=1,ng
c  90 cx=cx+cm(l,j)*tc(l)
c  80 write(*,'(1x,a,6f12.5)') 'compare',roots(i),y(i),cx
   10 continue
      deallocate(roots,tc,y,d)
      end
