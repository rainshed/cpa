      subroutine lsqftr(x,y,n,p,m)
c-----------------------------------------------------------------------
c     (m-1)th order least square fit using n data.
c     f(x)=p(1)*x**(m-1)+p(2)*x**(m-2)+...+p(m)
c     coded by H. Akai, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (nmx=10)
      real*8 x(n),y(n),a(nmx,nmx),q(nmx),p(m),t(nmx)
      if(m .gt. nmx) call errtrp(1,'lsqftr','n too large')
      do 10 i=1,m
      q(i)=0d0
      do 10 j=i,m
   10 a(i,j)=0d0
      do 20 i=1,n
      t(m)=1d0
      do 30 j=m-1,1,-1
   30 t(j)=t(j+1)*x(i)
      do 20 j=1,m
      q(j)=q(j)+y(i)*t(j)
      do 20 l=j,m
   20 a(j,l)=a(j,l)+t(l)*t(j)
      call rsymeq(a,nmx,q,p,m)
      return
      end
