      function lindx(j)
c-----------------------------------------------------------------------
c     Given a j-index, this returns corresponding l-index.
c     coded by H.Akai, 1983, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      lindx=sqrt(dble(max(0,j-1)))+1.00001d0
      return
      end
