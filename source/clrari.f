      subroutine clrari(ia,n)
c-----------------------------------------------------------------------
c     Clear integer array of length n.
c     coded by H.Akai, 1983, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ia(n)
      do 10 k=1,n
   10 ia(k)=0d0
      return
      end
