      subroutine biomix(v1,v2,a,meshr)
c-----------------------------------------------------------------------
c     Mix v1 and v2 with mixing parameter a.
c     coded by H.Akai, 1984, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension v1(meshr),v2(meshr)
      b=1d0-a
      do 10 k=1,meshr
   10 v2(k)=b*v1(k)+a*v2(k)
      return
      end
