      subroutine equari(ia,ib,n)
c-----------------------------------------------------------------------
c     Equates integer numbers.
c     coded by H.Akai, 1983, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer ia(n),ib(n)
      do 10 k=1,n
   10 ib(k)=ia(k)
      return
      end
