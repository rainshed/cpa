      subroutine equarr(a,b,n)
c-----------------------------------------------------------------------
c     Equates real numbers.
c     coded by H.Akai, 1983, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n),b(n)
      do 10 k=1,n
   10 b(k)=a(k)
      return
      end
