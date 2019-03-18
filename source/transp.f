      subroutine transp(a,jmx)
c-----------------------------------------------------------------------
c     transposes a given matrix a of rank jmx.
c     coded by H.Akai, July, 1992, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 a(jmx,jmx)
      do 10 i=1,jmx-1
      do 10 j=i+1,jmx
      swap=a(i,j)
      a(i,j)=a(j,i)
   10 a(j,i)=swap
      return
      end
