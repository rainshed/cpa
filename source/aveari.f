      subroutine aveari(ia,ib,n)
c-----------------------------------------------------------------------
c     give spin up/down averaged data
c     coded by H.Akai, 1993, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer ia(n),ib(n)
      do 10 k=1,n
      ia(k)=(ia(k)+ib(k))/2
   10 ib(k)=ia(k)
      return
      end
