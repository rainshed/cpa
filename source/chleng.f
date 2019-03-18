      subroutine chleng(a,ln)
c---------------------------------------------------------------------
c     Returns the length of the non blank part of the character a.
c     coded by H.Akai
c---------------------------------------------------------------------
      character a*(*)
      n=len(a)
      do 10 i=n,1,-1
      ln=i
      if(a(i:i) .ne. ' ') return
   10 continue
      ln=0
      return
      end
