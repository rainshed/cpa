      subroutine swparr(a,b,n)
c-----------------------------------------------------------------------
c     swap spin up/down data
c     coded by H.Akai, 1993, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 a(n),b(n)
      do 10 k=1,n
      swap=b(k)
      b(k)=a(k)
   10 a(k)=swap
      return
      end
