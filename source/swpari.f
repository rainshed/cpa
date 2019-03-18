      subroutine swpari(ia,ib,n)
c-----------------------------------------------------------------------
c     swap spin up/down data
c     coded by H.Akai, 1993, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer ia(n),ib(n)
      do 10 k=1,n
      iswap=ib(k)
      ib(k)=ia(k)
   10 ia(k)=iswap
      return
      end
