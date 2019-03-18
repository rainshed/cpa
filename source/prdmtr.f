      subroutine prdmtr(a,b,c,n)
c-----------------------------------------------------------------------
c     make the product of the n*n matrices a and b and put the
c     result into c.
c     coded by H.Akai, 14 Dec. 1996, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 a(n,n),b(n,n),c(n,n)
      do 10 i=1,n
      do 10 j=1,n
      c(j,i)=0d0
      do 10 k=1,n
   10 c(j,i)=c(j,i)+a(j,k)*b(k,i)
      end
