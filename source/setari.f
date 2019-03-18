      subroutine setari(ia,ib,n)
c-----------------------------------------------------------------------
c     Give a constant ib to an integer array ia.
c     coded by H.Akai, 1983, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer ia(n)
      do 10 k=1,n
   10 ia(k)=ib
      return
      end
