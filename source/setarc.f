      subroutine setarc(ca,cb,n)
c-----------------------------------------------------------------------
c     Give a constant cb to a complex*16 array ca.
c     coded by H.Akai, 1983, Juelich
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 ca(n),cb
      do 10 k=1,n
   10 ca(k)=cb
      return
      end
