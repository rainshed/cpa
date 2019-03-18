      subroutine prdmtc(a,b,c,n,m)
c-----------------------------------------------------------------------
c     make the product of the n*n complex matrices a and b and put the
c     result into c.
c     coded by H.Akai, 24 Nov. 1997, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 a(n,n),b(n,n),c(n,n)
      do 10 i=1,m
      do 10 j=1,m
      c(i,j)=(0d0,0d0)
      do 10 k=1,m
   10 c(i,j)=c(i,j)+a(i,k)*b(k,j)
      end
