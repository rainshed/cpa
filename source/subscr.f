      subroutine subscr(l,m,ml)
c-----------------------------------------------------------------------
c     Extract l amd m values from ml-index.
c     coded by H.Akai, 1983, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      l=sqrt(dble(max(0,ml-1)))+1d-5
      m=ml-l*(l+1)-1
      return
      end
