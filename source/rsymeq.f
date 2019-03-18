      subroutine rsymeq(a,m,b,x,n)
c-----------------------------------------------------------------------
c     Solves a set of real symmetric liner equations
c     using Gauss sweeping out with no pivoting
c     coded by H.Akai, 1975, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(m,n),b(n),x(n)
      do 10 k=1,n-1
      do 10 i=k+1,n
      c=a(k,i)/a(k,k)
      b(i)=b(i)-c*b(k)
      do 10 j=i,n
   10 a(i,j)=a(i,j)-c*a(k,j)
      x(n)=b(n)/a(n,n)
      do 20 i=n-1,1,-1
      x(i)=b(i)
      do 30 k=i+1,n
   30 x(i)=x(i)-a(i,k)*x(k)
   20 x(i)=x(i)/a(i,i)
      return
      end
