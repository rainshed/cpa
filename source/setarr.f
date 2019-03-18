      subroutine setarr(a,b,n)
c-----------------------------------------------------------------------
c     Give a constant b to a real*8 array a.
c     coded by H.Akai, 1983, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 a(n)
      do 10 k=1,n
   10 a(k)=b
      return
      end
