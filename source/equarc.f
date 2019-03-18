      subroutine equarc(a,b,n)
c-----------------------------------------------------------------------
c     Equate complex numbers.
c     coded by H.Akai, 1983, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n),b(n)
      complex*16 a,b
      do 10 k=1,n
   10 b(k)=a(k)
      return
      end
