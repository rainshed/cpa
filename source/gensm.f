      subroutine gensm(e,wt,f,mmxl,mxlcmp,ncmpx,kmx,sm,ng)
c-----------------------------------------------------------------------
c     Integrate f*t function, where t is the thebyceff polinomials,
c     and construct sm matrix which is readily used to calculate
c     charge density, etc.
c     coded by H.Akai, 1983, Juelich
c     modified by H. Akai, 1999, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 e(kmx),wt(ng,3,kmx),f(mmxl,ncmpx,kmx),a,b,c
      real*8 sm(ng,mmxl,ncmpx)
      integer mxlcmp(ncmpx)
      pi=4d0*atan(1d0)
      call clrarr(sm,ng*mmxl*ncmpx)
      do 10 k=1,kmx-2,2
      do 20 i=1,ncmpx
      do 20 l=1,mxlcmp(i)**2
      a=((f(l,i,k+1)-f(l,i,k))/(e(k+1)-e(k))-(f(l,i,k+2)-f(l,i,k+1))
     &   /(e(k+2)-e(k+1)))/(e(k)-e(k+2))
      b=(f(l,i,k+1)-f(l,i,k))/(e(k+1)-e(k))-a*(e(k+1)+e(k))
      c=f(l,i,k)-(a*e(k)+b)*e(k)
      do 20 j=1,ng
   20 sm(j,l,i)=sm(j,l,i)-dimag(a*(wt(j,3,k+2)-wt(j,3,k))
     &   +b*(wt(j,2,k+2)-wt(j,2,k))+c*(wt(j,1,k+2)-wt(j,1,k)))/pi
   10 continue
c     write(*,'(1x,1p,9e13.5)')(abs(f(9,1,k)),k=1,kmx)
c     write(*,'(1x,1p,6e13.5)')((sm(j,l,1),j=1,ng),l=1,9)
c  22 write(6,1100)k,e(k),(f(l,1,k),l=1,3),e(k+1),(f(l,1,k+1),l=1,3)
c    &            ,(sm(1,l,1),l=1,3)
c1100 format(/'                    k=',i3
c    &      ,2(/'  e=',2f10.5,'  f=',4(2f10.5,3x))
c    &       /'     sm=',1p,4e14.6)
      end
