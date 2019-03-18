      subroutine clrarc(a,n)
c-----------------------------------------------------------------------
c     Give zero value for complex array of length n.
c     coded by H.Akai, 1983, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n)
      complex*16 a
      do 10 k=1,n
   10 a(k)=(0d0,0d0)
      return
      end
