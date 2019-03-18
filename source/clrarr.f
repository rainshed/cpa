      subroutine clrarr(a,n)
c-----------------------------------------------------------------------
c     Clear real*8 arrey of length n.
c     coded by H.Akai, 1983, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n)
      do 10 k=1,n
   10 a(k)=0d0
      return
      end
